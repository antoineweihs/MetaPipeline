#' @title Save file functions
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @description functions that will save files as *.RData or *.csv
#' @param data (any) data set to be saved
#' @param type (string) either \code{csv} or \code{rdata}. File type \code{data} should be saved as
#' @param save_path (sting) full path to which the output should be saved to (including filename)
#' @param phenotype (string) phenotype used in the analysis
#' @param stratum (string) stratum used in the analysis
#' @param print_log (bool) TRUE: print to log file FALSE: won't
#' @param log_path (string) path to log file
#' @param verbose (bool) TRUE: Prints output to terminal FALSE: won't
#'
#' @importFrom utils write.csv
#'
#' @export
save_files <- function(data, type, save_path, phenotype, stratum, print_log, log_path, verbose)
{
  if(type == "csv")
  {
    #terminal and log file output
    text = paste0("Saving to (csv format): ", save_path)
    if(verbose) {writeLines(text)}
    if(print_log) {cat(c(text,""), file=log_path, append=TRUE, sep="\n")}

    #save csv data set
    utils::write.csv(data, file=save_path, quote=FALSE, row.names = FALSE)
  }
  else
  {
    #terminal and log file output
    text = paste0("Saving to (RData format): ", save_path)
    if(verbose) {writeLines(text)}
    if(print_log) {cat(c(text,""), file=log_path, append=TRUE, sep="\n")}

    #save as RData
    save(data, file=save_path)
  }
}
