#' @title post process for Z-based meta analysis
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @description post processing function for \code{Z} model which can draw
#'              draws a manhattan plot and/or annotates the output with an annotation file
#' @param result (data.frame) result from \code{\link{meta}} using Z
#' @param FDR (bool) TURE: FDR values are used instead of p-values (input needs to contain an FDR column) FALSE: p-values are used
#' @param significance_level (double) FDR/p-value significance level for cut off
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
#' @importFrom readr read_delim
#' @export
post_processZ <- function(result, FDR, significance_level=0.05, plot_manhattan=TRUE, output_path="./", annotate_result=TRUE,
                         annotation_filepath=NULL, phenotype=NULL, stratum=NULL, print_log=FALSE, log_path="./", verbose=TRUE)
{
  model="Z"
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
    relevant_id = result$Markername[which(as.numeric(as.character(result$Estimated_PValue)) < significance_level & result$analysis_error == FALSE)]
  }
  if(FDR == FALSE) {relevant_id = result$Markername[which(as.numeric(as.character(result$Estimated_PValue)) < significance_level & result$analysis_error == FALSE)]}

  #terminal and log file output
  text = paste0("Number of significant sites: ", length(relevant_id))
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  if(plot_manhattan & FDR)
  {
    text = "Plotting Manhattan plot with FDR values"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    temp_data = data.frame(Markername = result$Markername, Estimate_phenotype = rep(1, length(result$Estimated_PValue)), Estimated_PValue = result$Estimated_PValue)
    temp_data = merge(temp_data, Masterfile, by="Markername", all.x=TRUE, all.y=FALSE)
    test = double_manhattan(temp_data, chr="CHR", bp="MAPINFO", p="Estimated_PValue", markername="Markername", beta="Estimate_phenotype", FDRcorr=FDR, cutoff=significance_level,
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
    text = "Plotting Manhattan plot with P-values"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    temp_data = data.frame(Markername = result$Markername, Estimate_phenotype = rep(1, length(result$Estimated_PValue)), Estimated_PValue = result$Estimated_PValue)
    temp_data = merge(temp_data, Masterfile, by="Markername", all.x=TRUE, all.y=FALSE)
    test = double_manhattan(temp_data, chr="CHR", bp="MAPINFO", p="Estimated_PValue", markername="Markername", beta="Estimate_phenotype", FDRcorr=FDR, cutoff=significance_level,
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
