#' @title Load files for this package
#' @description This function imports the required files and modifies them in such a way, that they can be used by the package.
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_summary_path (string) path to the file that contains all the info of the cohorts. Has to include the following tab separated
#'                          columns: \itemize{
#'                          \item{cohort}{: Name of the cohort}
#'                          \item{phenotype}{: Name of the observed phenotype used in the file}
#'                          \item{stratum}{: Name of the stratum used in the file}
#'                          \item{array_type}{: type of array used by the cohort. Can be 450 or EPIC}
#'                          \item{probeID}{: Name of the column containing the probeIDs}
#'                          \item{beta}{: Name of the column containing the effect sizes}
#'                          \item{se}{: Name of the column containing the standard error values}
#'                          \item{p_val}{: Name of the column containing the p-values}
#'                          \item{size}{: Name of the column containing the sample sizes}
#'                          \item{separator}{: type of separator used in the file. Example: ","}
#'                          \item{filename}{: Full name of the file. Example: "example_name.csv"}
#'                          \item{file_path}{: Full path to the file. Example: "/home/user/file_location/"}
#'                          }
#' @param phenotype (string) phenotype that should be analysed (has to appear in the summary$phenotype column)
#' @param stratum (string) stratum that should be analysed (has to appear in the summary$stratum column)
#' @param cohort (string) if only certain cohorts of the summary file should be analysed (has to appear in the summary$cohort column)
#' @param FDR (bool) TRUE: Add a column to the output with Bonferroni corrected p-values FALSE: won't
#' @param annotation (bool) TRUE: Adds chromosome and MAPINFO column to output FALSE: won't
#' @param anno_file_path (string) path to annotation file (required if annotation is TRUE). NOTE:
#'                       - Has to include the following tab separated columns:
#'                         IlmnID (contains the CpG site IDs), CHR (contains the chromosome number), MAPINFO (contains the MAPINFO)
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log

#' @return data frame used by the rest of the package containing the combined data of the loaded cohorts

#' @importFrom data.table fread
#' @importFrom utils read.table
#' @import R.utils
#' @importFrom stats p.adjust

#' @export
load_files <- function(data_summary_path, phenotype, stratum, cohort=NULL, FDR=TRUE, annotation = FALSE,
                       anno_file_path= "./", verbose=TRUE, print_log=TRUE, log_path="./log.txt")
{
  ##output
  if(verbose) {writeLines("Loading data")}

  ## load summary of cohort data
  if(!file.exists(data_summary_path)) {stop(paste0("data summary file is not at ", data_summary_path, ". Please change the data_summary_path variable"))}
  data_summary = utils::read.table(data_summary_path, header = TRUE, sep = "\t", row.names= NULL)

  ##control step: check integrity of data_summary
  "%nin%" = Negate("%in%")
  data_summary_name = names(data_summary)
  column_names = c("phenotype", "stratum", "cohort", "file_path", "file_Name", "separator", "array_type")
  for (i in 1:length(column_names)) {if( column_names[i] %nin% data_summary_name) {stop(paste0(column_names[i], " column missing in data_summary"))}}

  ##control step: check if stratum and phenotype are present in data_summary
  if(stratum %nin% data_summary$stratum) {stop(paste0(stratum, " not appearing in the summary$stratum column"))}
  if(phenotype %nin% data_summary$phenotype) {stop(paste0(phenotype, " not appearing in the summary$phenotype column"))}
  if(!is.null(cohort) && !all(cohort %in% data_summary$cohort)) {stop("One of the cohort names not appearing in summary$cohort")}

  ## get data set info
  data_set_names = c()
  data_set_path = c()
  data_set_position = c()
  data_set_separator = c()
  data_set_cohort = c()
  data_set_array = c()

  for(i in 1:length(data_summary$file_Name))
  {
    if(!is.null(cohort))
    {
      if (data_summary$phenotype[i] == phenotype && data_summary$stratum[i] == stratum && data_summary$cohort[i] %in% cohort)
      {
        data_set_names = append(data_set_names, paste0(data_summary$cohort[i], "_", data_summary$phenotype[i], "_", data_summary$stratum[i]))
        data_set_path = append(data_set_path, paste0(data_summary$file_path[i], data_summary$file_Name[i]))
        data_set_separator = append(data_set_separator, as.character(data_summary$separator[i]))
        data_set_position = append(data_set_position, i)
        data_set_cohort = append(data_set_cohort, as.character(data_summary$cohort[i]))
        data_set_array = append(data_set_array, as.character(data_summary$array_type[i]))
      }
    }
    else
    {
      if (data_summary$phenotype[i] == phenotype && data_summary$stratum[i] == stratum)
      {
        data_set_names = append(data_set_names, paste0(data_summary$cohort[i], "_", data_summary$phenotype[i], "_", data_summary$stratum[i]))
        data_set_path = append(data_set_path, paste0(data_summary$file_path[i], data_summary$file_Name[i]))
        data_set_separator = append(data_set_separator, as.character(data_summary$separator[i]))
        data_set_position = append(data_set_position, i)
        data_set_cohort = append(data_set_cohort, as.character(data_summary$cohort[i]))
        data_set_array = append(data_set_array, as.character(data_summary$array_type[i]))
      }
    }
  }

  ## number of cohorts
  num_cohorts = length(data_set_names)

  ## load data set
  for (i in 1:num_cohorts)
  {
    text = paste0("   loading: ", data_set_names[i])
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    ## control step: check if file exists
    if(!file.exists(data_set_path[i])) {stop(paste0(data_set_names[i] , " could not be found at ", data_set_path[i]))}
    assign(data_set_names[i],  data.frame(data.table::fread(data_set_path[i], header= TRUE, sep= data_set_separator[i], verbose=FALSE)))
  }

  ##load an redefine the annotation file if required
  if(annotation)
  {
    anno_file = data.frame(data.table::fread(anno_file_path, sep="\t", verbose=FALSE))
    anno_file = data.frame(probeID=anno_file$IlmnID, CHR=anno_file$CHR, POS=anno_file$MAPINFO)
    anno_file$CHR = as.character(anno_file$CHR)
    anno_file$CHR[anno_file$CHR == "X"] = 23
    anno_file$CHR[anno_file$CHR == "Y"] = 24
  }


  ## rename and remove useless data
  if(verbose) {writeLines("   renaming and merging")}
  for(i in 1:num_cohorts)
  {
    temp_data = rename(get(data_set_names[i]), data_summary, data_set_position[i])
    if(FDR)
    {
      temp_data$FDR = stats::p.adjust(as.numeric(temp_data$P_VAL), "BH")
      temp_data = data.frame(probeID = temp_data$probeID, BETA=temp_data$BETA, SE=temp_data$SE, P_VAL=temp_data$P_VAL, Size=temp_data$Size,
                             Cohort=data_set_cohort[i], Array=data_set_array[i], FDR=temp_data$FDR, stringsAsFactors=FALSE)
    }
    else {temp_data = data.frame(probeID = temp_data$probeID, BETA=temp_data$BETA, SE=temp_data$SE, P_VAL=temp_data$P_VAL, Size=temp_data$Size,
                                 Cohort=data_set_cohort[i], Array=data_set_array[i], stringsAsFactors=FALSE)}
    if(i == 1) {combined_data = temp_data}
    else {combined_data = rbind(combined_data, temp_data)}
  }

  if(annotation)
  {
    combined_data = merge(combined_data, anno_file, by="probeID")
  }

  combined_data$BETA = as.numeric(as.character(combined_data$BETA))
  combined_data$SE = as.numeric(as.character(combined_data$SE))
  combined_data$P_VAL = as.numeric(as.character(combined_data$P_VAL))
  combined_data$Size = as.numeric(as.character(combined_data$Size))

  if(verbose){writeLines("") }
  return(combined_data)
}

rename <- function(data_set, data_summary, position)
{
  ## control step: summary integrity check
  "%nin%" = Negate("%in%")
  column_names = c("probeID", "beta", "se", "p_val", "size")
  data_column_names = c(as.character(data_summary$probeID[position]), as.character(data_summary$beta[position]), as.character(data_summary$se[position]),
                        as.character(data_summary$p_val[position]), as.character(data_summary$size[position]))

  for(i in 1:length(column_names))
  {
    if(column_names[i] %nin% names(data_summary)) {stop(paste0(column_names[i], " column missing in data_summary"))}
    if(data_column_names[i] %nin% names(data_set)) {stop(paste0("wrong ", column_names[i], " label in summary file for ", data_summary$cohort[position]))}
  }

  ## renaming
  names(data_set)[names(data_set) == data_summary$probeID[position]] = "probeID"
  names(data_set)[names(data_set) == data_summary$beta[position]] = "BETA"
  names(data_set)[names(data_set) == data_summary$se[position]] = "SE"
  names(data_set)[names(data_set) == data_summary$p_val[position]] = "P_VAL"
  names(data_set)[names(data_set) == data_summary$size[position]] = "Size"

  return(data_set)
}
