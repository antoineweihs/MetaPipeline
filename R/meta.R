#' @title Meta-analysis function
#' @description Function that performs the actual meta-analysis
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param input_data (data.frame) data set created by \code{\link{load_files}}
#' @param test_run (bool) FOR DEBUG TRUE: only the first \code{test_run_length} sites are analysed FALSE: everything is analysed
#' @param test_run_length FOR DEBUG (int) number of site to be run when \code{test_run} is TRUE
#' @param num_cores (int) number of cores that should be used (only works in Linux systems)
#' @param model (string) Estimator model for tau. Standard \code{FE}. Otherwise \code{BAYES} for a Bayesian approach, \code{REML} for a
#'              restricted maximum likelihood approach (not recommended for small number of cohorts)
#' @param cohort (vector string) order in which the cohorts should appear in the output. If null random order is used.
#' @param cpg_filter (vector string) a vector containing the CpG sites that should be analysed
#' @param FDR (bool) TRUE: p values will be corrected false discovery rate, FALSE: won't WARNING: If TRUE, only sites no errors will be returned
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log file
#' @param save_result (bool) TRUE: saves the result as *.RData FALSE: won't
#' @param save_dest (string) path to save spot
#' @param save_name (string) file name for the saved result
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#'
#' @return data.frame containing the results of the meta-analysis. If the model is FE or REML the output will have the following columns:
#'         \itemize{
#'         \item{\code{Markername}}{: CpG Site ID}
#'         \item{\code{Model}}{: Name of the model used (either FE or REML)}
#'         \item{\code{Estimate_phenotype}}{: estimated effect}
#'         \item{\code{SE_phenotype}}{: corresponding standard error}
#'         \item{\code{Pval_phenotype}}{: corresponding p-value}
#'         \item{\code{Zval_phenotype}}{: corresponding z-value}
#'         \item{\code{Heterogeneity}}{: estimated heterogeneity (NA if model is FE)}
#'         \item{\code{I2}}{: I squared parameter. Describes the percentage of variation across studies that is due to heterogeneity rather than by chance}
#'         \item{\code{number_of_cohorts}}{: number of cohorts available for the site}
#'         \item{\code{Cohort Names}}{: Multiple columns (one for each cohort) with + if the cohort effect is positive, - if it is negative, 0 if the effect is zero and ? if the no effect is present}
#'         \item{\code{sample_size}}{: total number of samples included in this site}
#'         \item{\code{analysis_error}}{: TRUE if a problem occurred, FALSE if none occurred}
#'         \item{\code{error_message}}{: contains error message if analysis_error is TRUE}
#'         \item{\code{FDR}}{: exists if FDR is TRUE and contains the FDR adjusted p-values}
#'         }
#'         If the model is BAYES the output will have the following columns
#'         \itemize{
#'         \item{\code{Markername}}{: CpG Site ID}
#'         \item{\code{Model}}{: Name of the model used (in this case BAYES)}
#'         \item{\code{MAP_effect}}{: maximum a-posteriori effect estimate}
#'         \item{\code{MAP_tau}}{: maximum a-posteriori heterogeneity estimate}
#'         \item{\code{bayes_factor_tau_eq_0}}{: Bayes factor for the H0 hypothesis tau = 0}
#'         \item{\code{bayes_factor_effect_eq_0}}{: Bayes factor for the H0 hypothesis effect = 0}
#'         \item{\code{minimum_bayes_factor_tau_eq_0}}{: minimum Bayes factor for the H0 hypothesis tau = 0}
#'         \item{\code{minimum_bayes_factor_effect_eq_0}}{: Bayes factor for the H0 hypothesis effect = 0}
#'         \item{\code{post_probability_tau_smaler_than_one}}{: posterior probability of the heterogeneity being smaller than 1, at which stage the heterogeneity can be considered to be extreme}
#'         \item{\code{post_probability_tau_bigger_than_one}}{: posterior probability of the heterogeneity being bigger than 1, at which stage the heterogeneity can be considered to be extreme}
#'         \item{\code{post_probability_effect_smaler_than_zero}}{: posterior probability of the effect being smaller than zero}
#'         \item{\code{post_probability_effect_bigger_than_zero}}{: posterior probability of the effect being bigger than zero}
#'         \item{\code{time_taken}}{: analysis time take for this site in seconds}
#'         \item{\code{number_of_cohorts}}{: number of cohorts available for the site}
#'         \item{\code{sample_size}}{: total number of samples included in this site}
#'         \item{\code{Cohort Names}}{: Multiple columns (one for each cohort) with + if the cohort effect is positive, - if it is negative, 0 if the effect is zero and ? if the no effect is present}
#'         \item{\code{analysis_error}}{: TRUE if a problem occurred, FALSE if none occurred}
#'         \item{\code{error_message}}{: contains error message if analysis_error is TRUE}
#'         }
#'
#' @importFrom metafor rma
#' @importFrom bayesmeta bayesmeta dhalfcauchy
#' @importFrom tictoc tic toc
#' @importFrom parallel mclapply
#' @importFrom dplyr arrange slice filter distinct
#' @importFrom pbmcapply progressBar pbmclapply
#' @importFrom stats sd
#'
#' @export
meta <- function(	input_data,
                  test_run=FALSE,
                  test_run_length=1000,
                  num_cores=1,
                  model="FE",
                  cohort = NULL,
                  cpg_filter=NULL,
                  FDR=TRUE,
                  print_log=FALSE,
                  log_path="log.txt",
                  save_result=FALSE,
                  save_dest="./",
                  save_name=NULL,
                  verbose=TRUE)
{
  ## control step: check if model is correct
  if (model != "REML" && model != "FE" && model != "BAYES")
  {
    stop(c("##########ERROR##########", "No known model chosen ",
           "Models are: REML (Restricted maximum likelihood), FE (Fixed effect) or BAYES", paste0("You chose: ", model) ))
  }

  ## control step: check if save_dest exists
  if(save_result && !dir.exists(save_dest)) {stop(paste0("output directory ", save_dest, " does not exist. Please change the output_path variable"))}

  ## control step: check integrity of input_data
  "%nin%" = Negate("%in%")
  data_names = names(input_data)
  column_names = c("probeID", "Cohort", "BETA", "SE", "P_VAL")
  for (i in 1:length(column_names)) {if( column_names[i] %nin% data_names) {stop(paste0(column_names[i], " column missing in input_data"))}}

  ## terminal and log file output and start timer
  if(verbose) {tictoc::tic()}
  text = c("Running meta-analysis", paste0("   start time: ",Sys.time()))
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  ## remove BETA and SE NAs from data
  input_data = input_data[!is.na(input_data$BETA),]
  input_data = input_data[!is.na(input_data$SE),]

  ## calculate effect mean and sd
  data_mean = mean(input_data$BETA)
  data_sd = stats::sd(input_data$BETA)

  ## get list of cohort names
  if(!is.null(cohort))
  {
    check = TRUE
    cohort_list = unique(input_data$Cohort)
    if (length(cohort_list) != length(cohort)) {check = FALSE}
    for (i in 1:length(cohort)) {if(cohort[i] %nin% cohort_list) {check=FALSE}}
    if(check==FALSE)
    {
      text = "Meta error: given cohort list does not match the cohorts found in the data. Using random cohort list."
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
    else {cohort_list=cohort}
  }
  if(is.null(cohort)) {cohort_list = unique(input_data$Cohort)}

  ## only considers sites in cpg_filter if given
  if(!is.null(cpg_filter)) {input_data = dplyr::filter(input_data, probeID %in% cpg_filter)}

  ## creating list that saves first and last position of each CpG site in the data set after having sorted it
  if(verbose) {writeLines("   Analysing data set")}
  input_data = dplyr::arrange(input_data, probeID)   #sort data
  input_data$probeID = as.character(input_data$probeID)
  start = rep(NA, length(input_data$probeID))   #create vector containing the start positions
  end = rep(NA, length(input_data$probeID))   #create vector containing the end positions
  i = pos = 1
  if(verbose) {progress_bar = pbmcapply::progressBar(max=length(input_data$probeID), style="ETA")}   #start progress bar if terminal output is wanted
  while(i <= length(input_data$probeID))
  {
    if(verbose) {progress_bar$up(i)}   #update progress bar if terminal output is wanted
    start[pos] = i  #save start position
    current_site = input_data$probeID[i]
    repeat{i = i+1; if(input_data$probeID[i] != current_site || is.na(input_data$probeID[i])) {break}}  #look for end position
    end[pos] = i-1 #save end position
    pos = pos+1
  }

  if(verbose) {close(progress_bar)}   #end progress bar if terminal output is wanted
  start = start[!is.na(start)]  #remove unwanted tail (vectors were defined to be too big)
  end = end[!is.na(end)]   #remove unwanted tail (vectors were defined to be too big)

  ## if test_run is TRUE: defining number of unique CpG sites
  if(test_run) {start = start[1:min(length(start), test_run_length)]; end = end[1:min(length(end), test_run_length)]}

  ## function performing the meta-analysis for a single site i when using the REML and FE model
  doOneMETAFOR <- function(i)
  {
    data = dplyr::slice(input_data, start[i]:end[i]) #extract required data from the combined data set
    sample_size = sum(data$Size)

    #if the CpG site only appears in 1 cohort, Na is returned for that site
    if(length(data$probeID) == 1)
    {
      #create cohort columns (shows beta direction of each cohort) + if the efect is positive, - if the effect is negative, 0 if the effect is zero and ? if the effect is missing
      included_cohorts = rep("?", length(cohort_list))
      if(data$BETA[1] > 0 && !is.na(data$BETA[1])) {included_cohorts[which(cohort_list == data$Cohort[1])] = "+"}
      else if (data$BETA[1] < 0 && !is.na(data$BETA[1])) {included_cohorts[which(cohort_list == data$Cohort[1])] = "-"}
      else if (data$BETA[1] == 0 && !is.na(data$BETA[1])) {included_cohorts[which(cohort_list == data$Cohort[1])] = "0"}
      #create output vector
      output = c(data$probeID[1], model, data$BETA[1], data$SE[1], data$P_VAL[1], NA, NA, NA, 1, included_cohorts, sample_size, TRUE, "Only one cohort")
    }

    #if the CpG site appears in more then 1 cohort an meta analysis is done on that site
    else
    {
      #create cohort columns (shows beta direction of each cohort) + if the efect is positive, - if the effect is negative, 0 if the effect is zero and ? if the effect is missing
      included_cohorts = rep("?", length(cohort_list))
      for(j in 1:length(data$probeID))
      {
        if(data$BETA[j] > 0 && !is.na(data$BETA[j])) {included_cohorts[which(cohort_list == data$Cohort[j])] = "+"}
        else if (data$BETA[j] < 0 && !is.na(data$BETA[j])) {included_cohorts[which(cohort_list == data$Cohort[j])] = "-"}
        else if (data$BETA[j] == 0 && !is.na(data$BETA[j])) {included_cohorts[which(cohort_list == data$Cohort[j])] = "0"}
      }

      # run meta-analysis using metafor package
      current_result = try(metafor::rma(yi=data$BETA, sei=data$SE, method=model), silent=TRUE)

      # check if error occured, if yes, return Na line with error code to output, if no return results
      if(inherits(current_result, "try-error")) {output = c(data$probeID[1], model, NA, NA, NA, NA, NA, NA, length(data$probeID), included_cohorts, sample_size, TRUE, current_result[1])}
      else {output = c(data$probeID[1], model, current_result$b, current_result$se, current_result$pval, current_result$zval, current_result$tau2, current_result$I2, length(data$probeID), included_cohorts, sample_size, FALSE, NA)}
    }
    return(output)
  }

  ## function performing the meta-analysis for a single site i when using the BAYES model
  doOneBAYES <- function(i)
  {
    #extract relevant data from combined data set
    data = dplyr::slice(input_data, start[i]:end[i])
    sample_size = sum(data$size)

    #create cohort columns (shows beta direction of each cohort) + if the efect is positive, - if the effect is negative, 0 if the effect is zero and ? if the effect is missing
    included_cohorts = rep("?", length(cohort_list))
    for(j in 1:length(data$probeID))
    {
      if(data$BETA[j] > 0 && !is.na(data$BETA[j])) {included_cohorts[which(cohort_list == data$Cohort[j])] = "+"}
      else if (data$BETA[j] < 0 && !is.na(data$BETA[j])) {included_cohorts[which(cohort_list == data$Cohort[j])] = "-"}
      else if (data$BETA[j] == 0 && !is.na(data$BETA[j])) {included_cohorts[which(cohort_list == data$Cohort[j])] = "0"}
    }

    #if the CpG site only appears in 1 cohort, Na is returned for that site
    if(length(data$probeID) == 1)
    {
      output = c(data$probeID[1], model, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 1, included_cohorts, sample_size, TRUE, "only one CpG site present")
    }

    #if the CpG site appears in more then 1 cohort an meta analysis is done on that site
    else
    {
      # run meta-analysis using the bayesmeta package (effect prior is standard normal distribution and tau prior is half-cauchy distribution (recommended by Roever))
      current_result = try(bayesmeta::bayesmeta(y=data$BETA, sigma = data$SE, tau.prior= function(t){bayesmeta::dhalfcauchy(t, scale=1)}, mu.prior = c("mean"=0, "sd"=1)),
                           silent=TRUE)

      # check if error occured, if yes, return Na line with error code to output, if no return results
      if(inherits(current_result, "try-error")) {output = c(data$probeID[1], model, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                                                            included_cohorts, sample_size, TRUE, current_result)}
      else
      {
        output = c(data$probeID[1], model, current_result$MAP[2,2], current_result$MAP[2,1],  current_result$bayesfactor[1,"tau=0"], current_result$bayesfactor[1,"mu=0"],
                   current_result$bayesfactor[2,"tau=0"], current_result$bayesfactor[2,"mu=0"],current_result$pposterior(tau=1),
                   1-current_result$pposterior(tau=1), current_result$pposterior(mu=0), 1 - current_result$pposterior(mu=0), current_result$init.time,
                   length(data$probeID), included_cohorts, sample_size, FALSE, NA)
      }
    }
    return(output)
  }


  ## control step: check system (mclapply only works on Linux) if not, the package will set the number of cores to 1
  if(Sys.info()["sysname"] != "Linux")
  {
    text = c("#####################", "Setting number of cores to 1 (Parallel runs only work on Linux systems)", "#####################")
    if(verbose) {writeLines(text)}
    if(print_log) {cat(c(text, ""), file=log_path, append=TRUE, sep="\n")}
    num_cores = 1
  }

  ## terminal and log file output
  text = c(paste0("   Model: ", model), paste0("   Number of Cores: ", num_cores), paste0("   Analysing ", length(start), " positions"),
           paste0("   Test run: ", test_run))
  text = append(text, c(paste0("   FDR: ", FDR), paste0("   Overall effect mean: ", data_mean), paste0("   Overall effect sd: ", data_sd)))
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  ## run analysis
  options(mc.cores = num_cores)
  if (model == "REML" || model == "FE")
  {
    #if verbose is true the package will display a progress bar, if not not terminal output will be generated
    if(verbose && num_cores > 1) {meta_result = do.call(rbind, pbmcapply::pbmclapply(1:length(start), doOneMETAFOR, mc.style="ETA"))}
    else {meta_result = do.call(rbind, parallel::mclapply(1:length(start), doOneMETAFOR))}

    #define column names of the output and reconfigure the output data.frame
    colnames(meta_result) = c("Markername", "Model", "Estimate_phenotype", "SE_phenotype", "Pval_phenotype", "Zval_phenotype", "Heterogeneity", "I2",
                              "number_of_cohorts", cohort_list, "sample_size", "analysis_error", "error_message")
    meta_result = data.frame(meta_result, stringsAsFactors = FALSE)

    # set heterogenity at NA if it FE model is used (as it is not considered by this model)
    if(model == "FE") {meta_result$Heterogeneity = NA}

    # calculated FDR adjusted p-values of the meta-analysis outcome (will remove any site where an error occured)
    if(FDR)
    {
      ## terminal and log file output
      text = paste0("   Applying FDR (removing ", length(which(meta_result$analysis_error == TRUE)), " positions because of errors)")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(c(text), file=log_path, append=TRUE, sep="\n")}

      meta_result = meta_result[meta_result$analysis_error == FALSE,]
      meta_result$FDR = p.adjust(as.numeric(meta_result$Pval_phenotype), "BH")
    }
  }
  else if(model == "BAYES")
  {
    #if verbose is true the package will display a progress bar, if not not terminal output will be generated
    if(verbose && num_cores > 1) {meta_result = do.call(rbind, pbmcapply::pbmclapply(1:length(start), doOneBAYES, mc.style="ETA"))}
    else {meta_result = do.call(rbind, parallel::mclapply(1:length(start), doOneBAYES))}

    #define column names of the output and reconfigure the output data.frame
    colnames(meta_result) = c("Markername", "Model" , "MAP_effect", "MAP_tau", "bayes_factor_tau_eq_0", "bayes_factor_effect_eq_0",  "minimum_bayes_factor_tau_eq_0",
                              "minimum_bayes_factor_effect_eq_0", "post_probability_tau_smaler_than_one", "post_probability_tau_bigger_than_one",
                              "post_probability_effect_smaler_than_zero", "post_probability_effect_bigger_than_zero", "time_taken", "number_of_cohorts",
                              cohort_list, "sample_size", "analysis_error", "error_message")
    meta_result = data.frame(meta_result, stringsAsFactors = FALSE)
  }

  ## saves the result if save_result = TRUE
  if(save_result)
  {
    if(verbose) {writeLines("   Saving result")}
    if(is.null(save_name)) {save_name = paste0("meta_analysis_result_",model,"_", format(Sys.Date(), "%d_%m_%Y"),".RData")}
    save_name = paste0(save_dest, save_name)
    save(meta_result, file=save_name)
  }

  ## terminal and log file output
  text = c("   Done with meta-analysis", paste0("   End time: ", Sys.time()))
  if(verbose) {writeLines(text)}
  if(print_log) {cat(c(text,""), file=log_path, append=TRUE, sep="\n")}
  if(verbose) {tictoc::toc()}

  return(meta_result)
}
