#' @title Apply BACON
#' @description This function uses the \code{\link[bacon]{bacon}} package to calculate the inflation and bias and modifies the input data if inflation is above 1
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param run (bool) TRUE: this function is applied, FALSE: this function won't be applied
#' @param print_log (bool) TRUE: print to log file FALSE: won't
#' @param log_path (string) path to log file
#' @param verbose (bool) TRUE: print log to command line FALSE: won't
#' @return \code{data_set} with BETA and SE values modified if BACON-inflation is above 1. NOTE: Sites where effect or standard error was NA have been removed
#'
#' @importFrom bacon bacon inflation bias es se
#'
#' @export
run_bacon <- function(data_set, run=TRUE, print_log=FALSE, log_path="log.txt", verbose=TRUE)
{
  if(run)
  {
    ##terminal and log file output
    text= c("Running BACON", "   Cohort parameters: ")
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    #extract list of cohort names from data set
    cohorts = unique(data_set$Cohort)
    for(i in 1:length(cohorts))
    {
      #extract data from the consdered cohort and remove na's
      temp_data = data_set[data_set$Cohort == cohorts[i],]
      temp_data = temp_data[!is.na(temp_data$BETA),]
      temp_data = temp_data[!is.na(temp_data$SE),]

      if(length(temp_data$BETA) == 0)
      {
        text = paste0("    ",cohorts[i], " skiped, no sites present")
        if(verbose) {writeLines(text)}
        if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
      }
      else
      {
        #run bacon function
        bc = bacon::bacon(NULL, temp_data$BETA, temp_data$SE)

        #check if inflation is above 1
        if(bacon::inflation(bc) == "NaN")
        {
          text = paste0("    ",cohorts[i], " skiped, bacon output is NaN")
          if(verbose) {writeLines(text)}
          if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
        }
        else
        {
          if(bacon::inflation(bc) > 1)
          {
            ##terminal and log file output
            text = paste0("    ",cohorts[i], "-> inflation: ", bacon::inflation(bc), ", bias: ", bacon::bias(bc)," (modified)")
            if(verbose) {writeLines(text)}
            if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

            #modify effects and standard errors in the combined data file
            temp_data$BETA = as.vector(bacon::es(bc))
            temp_data$SE = as.vector(bacon::se(bc))
            temp_data$P_VAL = as.vector(bacon::pval(bc))
            data_set = data_set[data_set$Cohort != cohorts[i],]
            data_set = rbind(data_set, temp_data)
          }
          else
          {
            #terminal and log file output
            text = paste0("    ",cohorts[i], "-> inflation: ", bacon::inflation(bc), ", bias: ", bacon::bias(bc))
            if(verbose) {writeLines(text)}
            if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
          }
        }
      }
    }

    ##terminal and log file output
    if(print_log) {cat("", file=log_path, append=TRUE, sep="\n")}
    if(verbose) {writeLines("")}
  }

  return(data_set)
}
