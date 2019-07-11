#' @title Meta-analysis pipeline
#' @description Main script that handles the meta analysis. The script is written in a modular fashion and multiple entry points have been added to
#'              expand the script where needed. Please note that before adding new scripts at an entry point, it is advisable to check the source code
#'              to make sure the script is compatible. The steps of the script are: load files -> pre-filter entry point -> run filter -> pre-bacon
#'              entry point -> run bacon -> pre-meta-analysis entry point -> run meta-analysis -> post-meta-analysis entry point -> post process
#'              -> final entry point -> save results
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param phenotype (string) phenotype that should be analysed (has to appear in the summary$Phenotype column)
#' @param gender (string) gender that should be analysed (has to appear in the summary$Gender column)
#' @param num_cores (int) number of cores that should be used for the meta analysis (only works on Linux)
#' @param model (string) Estimator model for tau. Standard \code{FE}. Otherwise \code{BAYES} for a Bayesian approach,
#'          		\code{REML} for a restricted maximum likelihood approach (not recommended for small number of cohorts (<10))
#' @param data_summary_path (string) path to the file that contains all the info of the cohorts. Has to include the following tab separated
#'                          columns: \itemize{
#'                          \item{\code{cohort}}{: Name of the cohort}
#'                          \item{\code{phenotype}}{: Name of the observed phenotype used in the file}
#'                          \item{\code{gender}}{: Name of the gender used in the file}
#'                          \item{\code{array_type}}{: type of array used by the cohort. Can be 450 or EPIC}
#'                          \item{\code{probeID}}{: Name of the column containing the probeIDs}
#'                          \item{\code{beta}}{: Name of the column containing the effect sizes}
#'                          \item{\code{se}}{: Name of the column containing the standard error values}
#'                          \item{\code{p_val}}{: Name of the column containing the p-values}
#'                          \item{\code{size}}{: Name of the column containing the sample sizes}
#'                          \item{\code{separator}}{: type of separator used in the file. Example: ","}
#'                          \item{\code{filename}}{: Full name of the file. Example: "example_name.csv"}
#'                          \item{\code{file_path}}{: Full path to the file. Example: "/home/user/file_location/"}
#'                          }
#' @param cohort (string vector) if only certain cohorts of the summary file should be analysed (has to appear in the summary$Cohort column)
#' @param cpg_filter (string vector) if only a certain CpG site shall be analysed
#' @param pre_bacon_step (bool) TRUE: runs script located at \code{pre_bacon_path} FALSE: won't
#' @param pre_bacon_path (string) path to script. \code{pre_bacon_step} has to be TRUE to run
#' @param apply_bacon (bool) PREPROCESS: TRUE: \code{\link[bacon]{bacon}} will be run on cohorts that meet certain criteria FALSE: \code{\link[bacon]{bacon}} won't be applied
#' @param pre_filter_step (bool) TRUE: runs script located at \code{pre_filter_path} FALSE: won't
#' @param pre_filter_path (string) path to script. \code{pre_filter_step} has to be TRUE to run
#' @param filter (bool) PREPROCESS: TRUE: \code{\link{run_filter}} will be applied FALSE: won't
#' @param filter_percentage (int) PREPROCESS: if \code{filter} is TRUE -> percentage of the maximal sample size that has to be present in the
#'                          sample site of each site
#' @param FDR (bool) TRUE: p values will be corrected for FDR, FALSE: won't
#' @param pre_meta_step (bool) TRUE: runs script located at \code{pre_meta_path} FALSE: won't
#' @param pre_meta_path (string) path to script. \code{post_meta_path} has to be TRUE to run
#' @param post_meta_step (bool) TRUE: runs script located at \code{pre_meta_path} FALSE: won't
#' @param post_meta_path (string) path to script. \code{post_meta_step} has to be TRUE to run
#' @param run_post_process (bool) TRUE: applies post_process function, FALSE: Won't
#'                         NOTE: Significant runtime increase if applied using the BAYES model and \code{run_posterior_Check} = TRUE
#' @param run_posterior_Check (bool) TRUE: will run a posterior predictive (see \link[bayesmeta]{pppvalue}) (NOTE: computationally intensive) FALSE: won't
#' @param post_sigLevel (double) POSTPROCESS: if \code{run_post_process} is TRUE and \code{model} is FE or REML-> FDR significance level for cut off
#' @param post_cutoff (double) POSTPROCESS: if \code{run_post_process} is TRUE and \code{model} is BAYES -> effect posterior probability cutoff
#' @param post_replicates (int) POSTPROCESS: if \code{run_post_process} is TRUE and \code{model} is BAYES -> number of replicates performed during the
#'                        posterior probability check (n >> 100 recommended)
#' @param plot_forest (bool) POSTPROCESS: TRUE: will create forest plot saved to output_path if run_post_process is TRUE FALSE: won't
#' @param annotate_result (bool) POSTPROCESS: TRUE: will merge output file with annoation file FALSE: won't
#' @param annotation_filepath (string) POSTPROCESS: file path to annotation file (annotation file has to contain a 'Markername' column, has to be a tab separated file,
#'                            comments have to be preceded by '#')
#' @param final_step (bool) TRUE: runs script located at \code{final_path} FALSE: won't
#' @param final_path (string) path to script. \code{final_step} has to be TRUE to run
#' @param print_log (bool) TRUE: log will be printed to output_path FALSE: no log will be printed
#' @param save_RData (bool) TRUE: result will be saved as .RData file to output_path, FALSE: result won't be saved
#' @param save_csvData (bool) TRUE: result will be saved as .csv file to output_path, FALSE: result won't be saved
#' @param output_path (string) path were the result should be saved
#' @param verbose (bool) TRUE: will print information to the terminal FALSE: won't
#'
#' @return data.frame containing the results of the pipeline (see \code{\link{meta}} and \code{\link{bayes_post_process}}
#' @export
run_meta <- function(	phenotype,
                      gender,
                      num_cores=1,
                      model="FE",
                      data_summary_path="/home/weihsa/data/Thyroid_Cohorts/Cohort_Summary/summary.txt",
                      cohort=NULL,
                      cpg_filter=NULL,
                      pre_bacon_step = FALSE, pre_bacon_path,
                      apply_bacon=TRUE,
                      pre_filter_step = FALSE, pre_filter_path,
                      filter=TRUE, filter_percentage=25,
                      FDR=TRUE,
                      pre_meta_step=FALSE, pre_meta_path,
                      post_meta_step=FALSE, post_meta_path,
                      run_post_process=FALSE, run_posterior_Check = TRUE, post_sigLevel = 0.05, post_cutoff=0.95, post_replicates=1000,
                      plot_forest=FALSE, annotate_result = FALSE, annotation_filepath=NULL,
                      final_step=FALSE, final_path,
                      print_log=TRUE,
                      save_RData=FALSE,
                      save_csvData=FALSE,
                      output_path="./",
                      verbose=TRUE)
{
  ## control step: check if everything is there
  if(!file.exists(data_summary_path)) {stop(paste0("data summary file is not at ", data_summary_path, ". Please change the data_summary_path variable"))}
  if((save_RData | save_csvData) & !dir.exists(output_path)) {stop(paste0("output directory ", output_path,
                                                                          " does not exist. Please change the output_path variable"))}

  ## create path variable for log file and write intro
  text = c("","######## Welcome to run_meta ########", "", paste0("Start Time: ", Sys.time()), "Global input variables: ", paste0("   Phenotype: ", phenotype),
           paste0("   gender: ", gender), paste0("   num_cores: ", num_cores), paste0("   model: ", model),
           paste0("   data_summary_path: ", data_summary_path), paste0("   cohort: ", paste(cohort, collapse=", ")),
           paste0("   cpg_filter: ", paste(cpg_filter, collapse=", ")), paste0("   pre_bacon_step: ", pre_bacon_step),
           if(pre_bacon_step) {paste0("   pre_bacon_path: ", pre_bacon_path)}, paste0("   apply_bacon: ", apply_bacon),
           paste0("   pre_filter_step: ", pre_filter_step), if(pre_filter_step) {paste0("   pre_filter_path: ", pre_filter_path)},
           paste0("   filter: ", filter), paste0("   FDR: ", FDR), paste0("   pre_meta_step: ", pre_meta_step),
           if(pre_meta_step) {paste0("   pre_meta_path: ", pre_meta_path)}, paste0("   post_meta_step: ", post_meta_step),
           if(post_meta_step) {paste0("   post_meta_path: ", post_meta_path)}, paste0("   run_post_process: ", run_post_process),
           if(run_post_process) {paste0("   plot_forest: ", plot_forest)}, if(run_post_process) {paste0("   annotate_result: ", annotate_result)},
           if(run_post_process) {paste0("   annotation_filepath: ", annotation_filepath)}, if(run_post_process) {paste0("   run_posterior_Check: ", run_posterior_Check)},
           if(run_post_process) {paste0("   post_sigLevel: ", post_sigLevel)}, if(run_post_process) {paste0("   post_cutoff: ", post_cutoff)},
           if(run_post_process) {paste0("   post_replicates: ", post_replicates)}, paste0("   final_step: ", final_step),
           if(final_step) {paste0("   final_path: ", final_path)}, paste0("   print_log: ", print_log), paste0("   save_RData: ", save_RData),
           paste0("   save_csvData: ", save_csvData), paste0("   output_path: ", output_path),"")

  # if terminal output is wanted prints intro to terminal
  if(verbose) {writeLines(text)}

  # if log output is wanted, prints log to terminal
  if(print_log)
  {
    log_path = paste0(output_path, model, "_", phenotype,"_",gender,"_",format(Sys.Date(), "%d_%m_%Y"),"_log.txt")
    cat(text, file=log_path, append=TRUE, sep="\n")
  }
  else {log_path=NULL}

  ## Meta-analysis pipeline:
  ## load files
  combined_data = load_files(data_summary_path=data_summary_path, phenotype=phenotype, gender=gender, cohort=cohort, FDR=FALSE, annotation=FALSE,
                             anno_file_path=NULL, verbose=verbose, print_log=print_log, log_path=log_path)

  ## pre process
  # entry point for the pre_filter step
  if(pre_filter_step)
  {
    if(file.exists(pre_filter_step)) {source(pre_filter_step)}
    else
    {
      text = paste0(pre_filter_step, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }
  # run the filter if filter variable is TRUE
  combined_data = run_filter(data_set=combined_data, percentage=filter_percentage, run=filter, print_log=print_log, log_path=log_path,
                             verbose=verbose)

  # entry point for the pre_bacon step
  if(pre_bacon_step)
  {
    if(file.exists(pre_bacon_path)) {source(pre_bacon_path)}
    else
    {
      text = paste0(pre_bacon_path, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }
  # apply bacon modification if apply_bacon variable is true
  combined_data = run_bacon(data_set=combined_data, run=apply_bacon, print_log=print_log, log_path=log_path, verbose=verbose)

  ## run meta-analyis
  #entry point for the pre meta-analysis step
  if(pre_meta_step)
  {
    if(file.exists(pre_meta_path)) {source(pre_meta_path)}
    else
    {
      text = paste0(pre_meta_path, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }
  # run meta-analysis
  meta_result = meta(input_data=combined_data, num_cores=num_cores, model=model, cohort = cohort, cpg_filter=cpg_filter, FDR=FDR,
                     print_log=print_log, log_path=log_path, verbose=verbose)

  # post meta-analysis entry point
  if(post_meta_step)
  {
    if(file.exists(post_meta_path)) {source(post_meta_path)}
    else
    {
      text = paste0(post_meta_path, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

  ## post process
  ##run FE or REML post process (creates forest plots of significant sites)
  if(run_post_process & (model=="FE" | model=="REML"))
  {
    ## terminal and log file output
    text = c(" ","Running post process")
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    meta_result = post_process(result=meta_result, combined_data=combined_data, model=model, FDR=FDR, significance_level=post_sigLevel, plot_forest=plot_forest, annotate_result=annotate_result,
                  annotation_filepath=annotation_filepath, output_path=output_path, phenotype=phenotype, gender=gender, print_log=print_log,  log_path=log_path, verbose=verbose)

    ## terminal and log file output
    text = c(" ")
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  }

  ##run BAYES post process (creates forest plots of sites above cut-off and/or applies a posterior predictive check)
  if(run_post_process & model=="BAYES")
  {
    meta_result = bayes_post_process(result=meta_result, combined_data=combined_data, run_posterior_Check=run_posterior_Check,
                                     cutoff=post_cutoff, replicates=post_replicates, plot_forest=plot_forest, output_path=output_path,
                                     phenotype=phenotype, gender=gender, print_log=print_log, log_path=log_path, verbose=verbose,
                                     num_cores=num_cores)

    ## terminal and log file output
    text = c(" ")
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  }

  #final entry point
  if(final_step)
  {
    if(file.exists(final_path)) {source(final_path)}
    else
    {
      text = paste0(final_path, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

  #save meta-analysis output as RData (if save_RData = TRUE)
  if(save_RData)
  {
    path=paste0(output_path, model, "_", phenotype,"_",gender,"_",format(Sys.Date(), "%d_%m_%Y"),".RData")
    save_files(meta_result, "rdata", path, phenotype, gender, print_log, log_path, verbose)
  }

  #save meta-analysis output as csv (if save_csvData = TRUE)
  if(save_csvData)
  {
    path=paste0(output_path, model, "_", phenotype,"_",gender,"_",format(Sys.Date(), "%d_%m_%Y"),".csv")
    save_files(meta_result, "csv", path, phenotype, gender, print_log, log_path, verbose)
  }

  return(meta_result)
}
