#' @title post process for FE or REML Random Effect meta analysis
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @description post processing function for \code{FE} and \code{REML} model which draws a forest plot for all site whose FDR adjusted p-value is above \code{significance_level}
#'
#' @param result (data.frame) result from \code{\link{meta}} using FE and including FDR
#' @param combined_data (data.frame) data used to create the above result from \code{\link{load_files}}
#' @param model (string) name of the model (can be either FE or REML)
#' @param FDR (bool) TURE: FDR values are used instead of p-values (input needs to contain an FDR column) FALSE: p-values are used
#' @param significance_level (double) FDR significance level for cut off
#' @param plot_forest (bool) TRUE: Draws a forest plot of the top hits FALSE: won't
#' @param plot_manhattan (bool) TRUE: will create a manhattan plot of the results FALSE: won't
#' @param output_path (string) output path for the forest plot
#' @param annotate_result (bool) TRUE: will merge output file with annoation file FALSE: won't
#' @param annotation_filepath (string) file path to annotation file (annotation file has to contain a 'Markername' column, has to be a tab separated file,
#'                            comments have to be preceded by '#')
#' @param phenotype (string) phenotype currently examined (for the forest plot name and label)
#' @param stratum (string) stratum currently examined (for the forest plot and label)
#' @param print_log (bool) TRUE: print to log file FALSE: won't
#' @param log_path (string) path to log file
#' @param verbose (bool) TRUE: Prints output to terminal FALSE: won't
#'
#' @return data.frame containing the original \code{result} data set if \code{annotate_result} = FALSE or the \code{result} data.frame merged with the annotation data.frame if
#'         \code{annotate_result} = TRUE.
#'
#' @importFrom metafor rma forest
#' @importFrom readr read_delim
#' @export
post_process <- function(result, combined_data, model, FDR, significance_level=0.05, plot_forest=TRUE, plot_manhattan=TRUE, output_path="./", annotate_result=TRUE,
                         annotation_filepath=NULL, phenotype=NULL, stratum=NULL, print_log=FALSE, log_path="./", verbose=TRUE)
{
  #terminal and log file output
  text = "Running post processing on results"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  ## get significant sites
  "%nin%" = Negate("%in%")
  if(FDR & ("FDR" %in% names(result))) {relevant_id = result$Markername[which(result$FDR < significance_level)]}
  if(FDR & ("FDR" %nin% names(result)))
  {
    text = "Post Process error: FDR = TRUE but no FDR column in data set. Using p-values."
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    relevant_id = result$Markername[which(as.numeric(as.character(result$Pval_phenotype)) < significance_level & result$analysis_error == FALSE)]
  }
  if(FDR == FALSE) {relevant_id = result$Markername[which(as.numeric(as.character(result$Pval_phenotype)) < significance_level & result$analysis_error == FALSE)]}

  #terminal and log file output
  text = paste0("Number of significant sites: ", length(relevant_id))
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  #check if significant sites have been found
  if(length(relevant_id) == 0)
  {
    text = "No significant site. Returning original data set"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  }
  else
  {
    for(i in 1:length(relevant_id))
    {
      #extract relevant sites from the cohort data set
      temp_data = combined_data[combined_data$probeID == relevant_id[i],]

      ## create forest plot
      if(plot_forest)
      {
        #rerun model
        temp = metafor::rma(yi=temp_data$BETA, sei=temp_data$SE, method=model, slab=temp_data$Cohort)
        #save forest plot
        png(paste0(output_path, phenotype, "_", stratum, "_", relevant_id[i], "_forestplot.png"))
        metafor::forest(temp, xlab=paste0("effect ", phenotype, " ", stratum, " ", relevant_id[i]))
        dev.off()
      }
    }
  }

  if(plot_manhattan & FDR)
  {
    text = "Plotting double Manhattan plot with FDR values"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    temp_data = data.frame(Markername = result$Markername, Estimate_phenotype = result$Estimate_phenotype, Pval_phenotype = result$Pval_phenotype)
    temp_data = merge(temp_data, Masterfile, by="Markername", all.x=TRUE, all.y=FALSE)
    test = double_manhattan(temp_data, chr="CHR", bp="MAPINFO", p="Pval_phenotype", markername="Markername", beta="Estimate_phenotype", FDRcorr=FDR, cutoff=significance_level,
                            strict_cutoff=NULL, title=paste0(phenotype, " ", stratum, " Manhattan plot"), logP=T, marktop=T, save_plot=T, save_path=output_path)
    if(test == 1)
    {
      text = "Post Process error: result vector is empty. Skipping manhatan plot"
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }
  if(plot_manhattan & FDR == FALSE)
  {
    text = "Plotting double Manhattan plot with P-values"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    temp_data = data.frame(Markername = result$Markername, Estimate_phenotype = result$Estimate_phenotype, Pval_phenotype = result$Pval_phenotype)
    temp_data = merge(temp_data, Masterfile, by="Markername", all.x=TRUE, all.y=FALSE)
    test = double_manhattan(temp_data, chr="CHR", bp="MAPINFO", p="Pval_phenotype", markername="Markername", beta="Estimate_phenotype", FDRcorr=FDR, cutoff=significance_level,
                            strict_cutoff=NULL, title=paste0(phenotype, " ", stratum, " Manhattan plot"), logP=T, marktop=T, save_plot=T, save_path=output_path)
    if(test == 1)
    {
      text = "Post Process error: result vector is empty. Skipping manhatan plot"
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
  }

  if(annotate_result)
  {
    text = "Annotating results"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    if(is.null(annotation_filepath))
    {
      text = "Post Process error: No annotation file given, skipping annotation step"
      if(verbose) {writeLines(text)}
      if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    }
    else
    {
      annotation_file <- readr::read_delim(annotation_filepath, "\t", escape_double = FALSE, comment = "#", trim_ws = TRUE, col_types = readr::cols(), progress = FALSE)
      result = merge(result, annotation_file, by="Markername", all.x=TRUE, all.y=FALSE)
    }
  }

  return(result)
}
