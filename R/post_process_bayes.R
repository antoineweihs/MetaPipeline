#' @title Bayes meta-analysis post process
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @description post processing function, that, if(\code{plot_forest} is TRUE plots forest plots of the sites whose effect posterior is above \code{cutoff} and,
#'              if \code{run_posterior_check} is TRUE, runs a posterior predictive check on the results. NOTE: Beware of significant runtime
#' @param result (data.frame) result from \code{\link{meta}} using the Bayesian model
#' @param combined_data (data.frame) data used to create the above result from \code{\link{load_files}}
#' @param run_posterior_Check (bool) TRUE: will run a posterior predictive (see \link[bayesmeta]{pppvalue}) (NOTE: computationally intensive) FALSE: won't
#' @param cutoff (double) effect posterior probability cut-off
#' @param replicates (int) number of replications for the posterior probability check (n >> 100 recommended)
#' @param plot_forest (bool) TRUE: Draws a forest plot of the top hits FALSE: won't
#' @param output_path (string) output path for the forest plot
#' @param annotate_result (bool) TRUE: will merge output file with annoation file FALSE: won't
#' @param annotation_filepath (string) file path to annotation file (annotation file has to contain a 'Markername' column, has to be a tab separated file,
#'                            comments have to be preceded by '#')
#' @param phenotype (string) phenotype used in the analysis
#' @param gender (string) gender used in the analysis
#' @param print_log (bool) TRUE: print to log file FALSE: won't
#' @param log_path (string) path to log file
#' @param verbose (bool) TRUE: Prints output to terminal FALSE: won't
#' @param num_cores (int) number of parallel runs (when \code{NULL}, the number of cores will be calculated automatically)
#'
#' @return \code{result} data frame containing all analysed sites (see \code{\link{meta}}) with the following extra columns
#'         \itemize{
#'         \item{\code{pp.applied}}{: TRUE if posterior check has been applied, FALSE else}
#'         \item{\code{pp.pValue}}{: p-value posterior predictive check (NA if effect posterior probability below \code{cutoff}) }
#'         \item{\code{pp.replicates}}{: number of replicates performed (equals the \code{replicates} variable)}
#'         \item{\code{pp.computation_time}}{: time required for the computation}
#'         }
#'
#' @importFrom bayesmeta bayesmeta forestplot.bayesmeta pppvalue dhalfcauchy
#' @importFrom grDevices dev.off png
#'
#' @export
bayes_post_process <- function(result, combined_data, run_posterior_Check = TRUE, cutoff=0.95, replicates=1000, plot_forest=TRUE, output_path="./",
                               annotate_result=TRUE, annotation_filepath=NULL,
                               phenotype=NULL, gender=NULL, print_log=FALSE, log_path="./", verbose=TRUE, num_cores=NULL)
{
  #terminal and log file output
  text = "Running Bayes post process"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  #identify sites with posterior effect probabilities above the cut-off
  significant_hits = which(result$post_probability_effect_bigger_than_zero >= cutoff | result$post_probability_effect_smaler_than_zero >= cutoff)

  #check if sites above cut-off exist, if no end function
  if(length(significant_hits)==0)
  {
    text = "No significant site. Returning original data set"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    return(result)
  }

  #terminal and log file output
  text = paste0("Found ", length(significant_hits), " sites with a probability above ", cutoff)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  #if posterior check is run add relevant columns to output
  if(run_posterior_Check)
  {
    result$pp.applied = FALSE
    result$pp.pValue = NA
    result$pp.replicates = NA
    result$pp.computation_time = NA
  }

  #go through sites abov cut-off
  for(i in 1:length(significant_hits))
  {
    #terminal and log file output
    text = paste0("Analysing site: ", result$Markername[significant_hits[i]])
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    #extract the relevant data from the cohort data set and rerun the meta-analysis
    data = combined_data[combined_data$probeID == result$Markername[significant_hits[i]],]
    current_result = try(bayesmeta::bayesmeta(y=data$BETA, sigma = data$SE, labels=data$Cohort, tau.prior= function(t){bayesmeta::dhalfcauchy(t, scale=1)},
                                              mu.prior = c("mean"=0, "sd"=1)), silent=TRUE) #half Cauchy prior recommended by Roever

    # check if the meta-analysis encountered problems if yes, skip to next site
    if(inherits(current_result, "try-error"))
    {
      text = paste0("Bayesmeta encountered an error at site ", result$Markername[significant_hits[i]], ". Error message: ", current_result)
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
    else
    {
      # if forest plot is wanted, creates a forest plot from the bayesmeta output
      if(plot_forest)
      {
        output_string = paste0(output_path, phenotype, "_", gender, "_", result$Markername[significant_hits[i]], "_forestplot.png")

        #terminal and log file output
        text = paste0("creating forest plot and saving it as: ", output_string)
        if(verbose) {writeLines(text)}
        if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

        #save output to file
        grDevices::png(output_string)
        bayesmeta::forestplot.bayesmeta(current_result, digits=8)
        grDevices::dev.off()
      }

      # if posterior predictive check is wanted, runs a one-sided check (significant run time increase)
      if(run_posterior_Check)
      {
        #checks which side the alternative hypothesis lies on
        myalternative = "greater"
        if(sign(as.numeric(as.character(result$MAP_effect[significant_hits[i]]))) == -1) {myalternative = "less"}

        #defines number of cores if undefined
        if(!is.null(num_cores)){num_cores=1}

        #runs posterior predictive check
        ppc_result = try(bayesmeta::pppvalue(current_result, parameter="mu", value=0, alternative=myalternative, statistic="cdf", n=replicates,
                                             parallel=num_cores, quietly=!verbose), silent=TRUE)

        #checks if error occured (sometimes due to the model running incorrectly), if yes tries to run the check again
        if(inherits(ppc_result, "try-error"))
        {
          #terminal and log file output
          text = paste0("Site ", result$Markername[significant_hits[i]], " encountered an error during the posterior predictive check. Trying again. Error message: ",
                        ppc_result)
          if(verbose) {writeLines(text)}
          if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

          #runs posterior predictive check again
          ppc_result = try(bayesmeta::pppvalue(current_result, parameter="mu", value=0, alternative=myalternative, statistic="cdf", n=replicates,
                                               if(!is.null(num_cores)){parallel=num_cores}, quietly=!verbose), silent=TRUE)
        }

        #checks if the check ran into problems again
        if(inherits(ppc_result, "try-error"))
        {
          #terminal and log file output
          text = paste0("Site ", result$Markername[significant_hits[i]], " again encountered an error during the posterior predictive check. Returning NA.
                        Error message: ", ppc_result)
          if(verbose) {writeLines(text)}
          if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

          # creates relevant output
          result$pp.applied[significant_hits[i]] = TRUE
          result$pp.pValue[significant_hits[i]] = NA
          result$pp.replicates[significant_hits[i]] = NA
          result$pp.computation_time[significant_hits[i]] = NA
        }
        #if no problems are present, saves the results to the output
        else
        {
          result$pp.applied[significant_hits[i]] = TRUE
          result$pp.pValue[significant_hits[i]] = ppc_result$p.value
          result$pp.replicates[significant_hits[i]] = as.numeric(ppc_result$parameter)
          result$pp.computation_time[significant_hits[i]] = ppc_result$computation.time
        }
      }
    }
  }

  if(annotate_result)
  {
    text = "Annotating output file"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    if(is.null(annotation_filepath))
    {
      text = "No annotation file given, skipping annotation step"
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
    else
    {
      text = "Annotating result file"
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
      annotation_file <- readr::read_delim(annotation_filepath, "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE)
      result = merge(result, annotation_file, by="Markername")
    }
  }

  return(result)
}
