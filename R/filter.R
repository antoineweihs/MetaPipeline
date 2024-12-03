#' @title Filter data
#' @description  This function filters the data, removing all sites that have a sample size below the threshold value
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param percentage (int) percentage of the maximal sample size that has to be present in the sample site of each site
#' @param run (bool) TRUE: this function is applied, FALSE: this function won't be applied
#' @param print_log (bool) TRUE: print to log file FALSE: won't
#' @param log_path (string) path to log file
#' @param verbose (bool) TRUE: print log to command line FALSE: won't
#'
#' @return filtered \code{data_set} NOTE: Sites with NA in effect or standard error column have been removed if run = T otherwise the initial dataset.
#'
#' @export
run_filter <- function(data_set, percentage=25, run=TRUE, print_log=FALSE, log_path="./log.txt", verbose=TRUE)
{
  if(run)
  {
    ##terminal and log file output
    text = "Applying size filter"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    ##extract cohort names from data set
    cohorts = unique(data_set$Cohort)
    for(i in 1:length(cohorts))
    {
      #extract cohort data from combined file and remove NA's
      temp_data = data_set[data_set$Cohort == cohorts[i],]
      temp_data = temp_data[!is.na(temp_data$BETA),]
      temp_data = temp_data[!is.na(temp_data$SE),]
      data_size_before = length(temp_data$Cohort)
      #calculate cut-off
      cutoff = (max(temp_data$Size, na.rm=T)/100)*percentage
      #remove sites where sample size is below the cut-off
      temp_data = temp_data[which(temp_data$Size > cutoff),]
      #add filtered data to the combined data
      data_set = data_set[data_set$Cohort != cohorts[i],]
      data_set = rbind(data_set, temp_data)

      #terminal and log file output
      text = paste0("   ", cohorts[i], ": Removed ", data_size_before - length(temp_data$Cohort), " Sites (cut-off: ", cutoff, ")")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }

    #terminal and log file output
    text=""
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  }
  return(data_set)
}
