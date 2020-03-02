#' @title Pre meta-analysis QC pipeline
#' @description This function combines multiple functions into one to perform an pre meta-analysis quality control. The script is written in a modular
#'              fashion and multiple entry points have been added to expand the script where needed. Please note that before adding new scripts at an
#'              entry point, it is advisable to check the source code it is advisable to check the source code to make sure the script is compatible. The
#'              steps of the script are: load_files -> pre_filter entry point -> run_filter -> pre_plotting entry point -> plotting -> plotting chromosomes
#'              -> pre_printing entry point -> printing summary -> printing flags -> final entry point
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_summary_path (string) path to the file that contains all the info of the cohorts. Has to include the following tab separated
#'                          columns: \itemize{
#'                          \item{\code{cohort}}{: Name of the cohort}
#'                          \item{\code{phenotype}}{: Name of the observed phenotype used in the file}
#'                          \item{\code{stratum}}{ Name of the stratum used in the file}
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
#' @param save_path (string) place where the outputs are saved
#' @param phenotype (string) phenotype that should be analysed (has to appear in the summary$Phenotype column)
#' @param stratum (string) stratum that should be analysed (has to appear in the summary$stratum column)
#' @param cohort (string) if only certain cohorts of the summary file should be analysed (has to appear in the summary$cohort column)
#' @param pre_filter_step (bool) TRUE: runs script located at \code{pre_filter_script} FALSE: won't
#' @param pre_filter_script (string) path to script. \code{pre_filter_step} has to be TRUE to run
#' @param run_filter (bool) TRUE: runs size filter, removing sites with a sample size of less then \code{filter_percentage} of the maximum sample size
#'                   FALSE: won't
#' @param filter_percentage (int) percentage used to filter
#' @param pre_plotting_step (bool) TRUE: runs script located at \code{pre_plotting_script} FALSE: won't
#' @param pre_plotting_script (string) path to script. \code{pre_plotting_step} has to be TRUE to run
#' @param plotting (bool) TRUE: plots will be generated FALSE: plots won't be generated
#' @param plotting_chr (bool) TRUE: chromosome plots will be generated FALSE: won't
#' @param same_scale (bool) TRUE: PVAL, BETA and SE separate plots are plotted using the same scale on the axes. FALSE: the separate plots are drawn on
#'                   different scale axes
#' @param plot_trim (bool) TRUE: will trim the density plots at local min and max FALSE: will plot density curve from global min to global max
#' @param pre_printing_step (bool) TRUE: runs script located at \code{pre_printing_script} FALSE: won't
#' @param pre_printing_script (string) path to script. \code{pre_printing_step} has to be TRUE to run
#' @param print_summary (bool) TRUE: file containing cohort summaries will be generated False: won't
#' @param print_flags (bool) TRUE: file containing flags will be generated FALSE: won't
#' @param flags_num_NA (int) parameter for \code{print_flags()}
#' @param flags_min_beta (int) parameter for \code{print_flags()}
#' @param flags_max_beta (int) parameter for \code{print_flags()}
#' @param flags_max_se (int) parameter for \code{print_flags()}
#' @param flags_min_inflation (int) parameter for \code{print_flags()}
#' @param flags_max_inflation (int) parameter for \code{print_flags()}
#' @param flags_decimal_places (int) parameter for \code{print_flags()}
#' @param final_step (bool) TRUE: runs script located at \code{final_script} FALSE: won't
#' @param final_script (string) path to script. \code{final_step} has to be TRUE to run
#' @param print_log (bool) TRUE: print to log file FALSE: won't
#' @param log_path (string) path to log file
#' @param verbose (bool) TRUE: print log to command line FALSE: won't
#'
#'@export
preMeta_QC<-function( data_summary_path,
                      save_path = "./", phenotype, stratum, cohort = NULL,
                      pre_filter_step = FALSE, pre_filter_script,                                                 #pre-filter steps(optional)
                      run_filter=TRUE, filter_percentage=25,                                                      #filter variables
                      pre_plotting_step = FALSE, pre_plotting_script,                                             #pre-plotting steps (optional)
                      plotting=TRUE, plotting_chr=TRUE,  same_scale=TRUE, plot_trim=FALSE,                        #plotting variables
                      pre_printing_step = FALSE, pre_printing_script,                                             #pre-printing steps (optional)
                      print_summary=TRUE,                                                                         #summary printing variables
                      print_flags=TRUE, flags_num_NA=1000, flags_min_beta=-2, flags_max_beta=2, flags_max_se=4,   #flag variables
                      flags_min_inflation=0.9, flags_max_inflation=1.1, flags_decimal_places=4,
                      final_step = FALSE, final_script,                                                           #final steps (optional)
                      print_log=FALSE, log_path="./log.txt", verbose=TRUE)                                        #log and terminal output variables

{
  # add a / at the end of the save_path string if missing
  if(substr(save_path, nchar(save_path), nchar(save_path)) != "/") {save_path = paste0(save_path, "/")}

  #control steps
  if(!file.exists(data_summary_path)) {stop(paste0("data summary file is not at ", data_summary_path, ". Please change the data_summary_path variable"))}
  if(!dir.exists(save_path)) {stop(paste0("output directory ", save_path, " does not exist. Please change the save_path variable"))}

  ##print output run parameters
  text = c("",paste0("### Analysing ", phenotype, " ", stratum, " ###"),"", "Input parameters:", paste0("   data_summary_path = ", data_summary_path),
           paste0("   save_path = ", save_path), paste0("   phenotype = ", phenotype), paste0("   stratum = ", stratum),
           if(!is.null(cohort)){paste0("   cohort = ", paste(cohort, collapse=", "))},
           paste0("   pre_filter_step = ", pre_filter_step), if(pre_filter_step) {paste0("   pre_filter_script = ", pre_filter_script)},
           paste0("   run_filter = ", run_filter), if(run_filter){paste0("   filter_percentage = ", filter_percentage)},
           paste0("   pre_plotting_step = ", pre_plotting_step), if(pre_plotting_step) {paste0("   pre_plotting_script = ", pre_plotting_script)},
           paste0("   plotting = ", plotting), paste0("   plotting_chr = ", plotting_chr), if(plotting) {paste0("   same_scale = ", same_scale)},
           if(plotting){paste0("   plot_trim = ", plot_trim)}, paste0("   pre_printing_step = ", pre_printing_step),
           if(pre_printing_step) {paste0("   pre_printing_script = ", pre_printing_script)}, paste0("   print_summary = ", print_summary),
           paste0("   print_flags = ", print_flags))
  if(print_flags) {text = c(text, paste0("   flag_num_NA = ", flags_num_NA), paste0("   flag_min_beta = ", flags_min_beta),
                            paste0("   flags_max_beta = ", flags_max_beta), paste0("   flags_max_se = ", flags_max_se),
                            paste0("   flags_min_inflation = ", flags_min_inflation), paste0("   flags_max_inflation = ", flags_max_inflation),
                            paste0("   flags_decimal_places = ", flags_decimal_places))}
  text = c(text, paste0("   final_step = ", final_step), if(final_step) {paste0("   final_script = ", final_script)},
           paste0("   print_log = ", print_log), if(print_log) {paste0("   log_path = ", log_path)}, paste0("   verbose = ", verbose))

  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  ## load_files
  combined_data <- load_files(data_summary_path=data_summary_path, phenotype=phenotype, stratum=stratum,
                              cohort=cohort, FDR=FALSE, annotation=TRUE,
                              verbose=verbose, print_log=print_log, log_path=log_path)

  ## optional pre_filter step
  if(pre_filter_step)
  {
    if(file.exists(pre_filter_script)) {source(pre_filter_script)}
    else
    {
      text = paste0(pre_filter_script, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

  # run filter
  if(run_filter) {combined_data <- run_filter(data_set=combined_data, percentage=filter_percentage, run=TRUE,
                                          log_path=log_path, verbose=verbose, print_log=print_log)}

  ## optional pre_plottig step
  if(pre_plotting_step)
  {
    if(file.exists(pre_plotting_script)) {source(pre_plotting_script)}
    else
    {
      text = paste0(pre_plotting_script, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

  ## create plots
  if(plotting)
  {
    beta_plot(data_set=combined_data,same_scale=same_scale, save_path=save_path, stratum=stratum, phenotype=phenotype,
              trim=plot_trim, verbose=verbose, print_log=print_log, log_path=log_path)
    pval_plot(data_set=combined_data,same_scale=same_scale, save_path=save_path, stratum=stratum, phenotype=phenotype,
              trim=plot_trim, verbose=verbose, print_log=print_log, log_path=log_path)
    se_plot(data_set=combined_data,same_scale=same_scale, save_path=save_path, stratum=stratum, phenotype=phenotype,
            trim=plot_trim, verbose=verbose, print_log=print_log, log_path=log_path)
    se_vs_size_plot(data_set=combined_data, save_path=save_path, stratum=stratum, phenotype=phenotype, verbose=verbose,
                    print_log=print_log, log_path=log_path)
  }
  if(plotting_chr)
  {
    chromosome_plot(data_set=combined_data,save_path=save_path, stratum=stratum, phenotype=phenotype, verbose=verbose,
                    print_log=print_log, log_path=log_path)
  }

  ## optional pre_printing step
  if(pre_printing_step)
  {
    if(file.exists(pre_printing_script)) {source(pre_printing_script)}
    else
    {
      text = paste0(pre_printing_script, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

  ## print summaries and flags
  if(print_summary) {print_summaries(data_set=combined_data, save_path=save_path, stratum=stratum, phenotype=phenotype,
                                     verbose=verbose, print_log=print_log, log_path=log_path)}
  if(print_flags) {print_flags(data_set=combined_data, num_NA=flags_num_NA, min_beta=flags_min_beta, max_beta=flags_max_beta,
                               max_se=flags_max_se, min_inflation=flags_min_inflation, max_inflation=flags_max_inflation,
                               decimal_places=flags_decimal_places, save_path=save_path, stratum=stratum, phenotype=phenotype,
                               verbose=verbose, print_log=print_log, log_path=log_path)}


  ## optional final step
  if(final_step)
  {
    if(file.exists(final_script)) {source(final_script)}
    else
    {
      text = paste0(final_script, " does not exists. Skipping this step")
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

}
