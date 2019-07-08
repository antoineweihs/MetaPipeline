#' @title Print Summaries
#' @description  This function prints summary statistics to an output file. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @export
print_summaries <- function(data_set,save_path, gender, phenotype, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  #sets number of decimal places to 12 digits
  options(digits=12)

  #define save path
  full_path = paste0(save_path,phenotype,"_",gender,"_cohort_stats.txt")

  #terminal and log file output
  text = paste0("printing: summaries to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  #create summary file
  sink(file= full_path, append= TRUE)
  cat("Summaries for ", phenotype, " ", gender, "\n", "\n")

  #extract cohort names from data set
  cohort_names = as.character(unique(data_set$Cohort))

  #create summaries for each cohort of the relevant columns in the data set
  for (i in 1:length(cohort_names))
  {
    temp_data = data_set[data_set$Cohort == cohort_names[i],]
    temp_data = data.frame(BETA=temp_data$BETA, SE=temp_data$SE, P_VAL=temp_data$P_VAL, SIZE=temp_data$Size)
    cat(cohort_names[i], "\n")
    print(summary(temp_data), quote=FALSE)
    se = temp_data$SE[!is.na(temp_data$BETA)]
    beta = temp_data$BETA[!is.na(temp_data$BETA)]
    beta = beta[!is.na(se)]
    se = se[!is.na(se)]
    cat("BACON inflation: ", bacon::inflation(bacon::bacon(NULL, beta, se)), " BACON bias: ", bacon::bias(bacon::bacon(NULL, beta, se)), "\n", "\n")
    #chisq = qchisq(1-temp_data$P_VAL, 1)
    #chisq = median(chisq, na.rm= TRUE)/qchisq(0.5,1)
    #cat("Lambda_GC: ", chisq,"\n","\n")
  }
  sink(file= NULL)
}


#' @title Print flags
#' @description  This function prints flags to an output file. Used for the pre-meta quality control analysis.
#'               At the moment these flags are:
#'               \itemize{
#'               \item{Different number of NAs in BETA, Pval and SE}
#'               \item{Number of NAs above \code{num_NA}}
#'               \item{min(BETA) below \code{min_beta} or max(BETA) above \code{max_beta}}
#'               \item{max(SE) above \code{max_se}}
#'               \item{more then 1 SE or BETA value above 75th + 1.5*IQR or below 25th percentile - 1.5*IQR}
#'               \item{\code{\link[bacon]{bacon}} inflation below \code{min_inflation} or above \code{max_inflation}}
#'               \item{number of cites of one chromosome = 0}
#'               \item{negative standard errors or p values}
#'               \item{average number of decimal places of BETA, Pval and SE smaller then \code{decimal_places}}
#'               }
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param num_NA (int) threshold for the number of NAs
#' @param min_beta (int) lower threshold for BETA values
#' @param max_beta (int) upper threshold for BETA values
#' @param max_se (int) upper threshold for SE values
#' @param min_inflation (int) lower threshold values for \code{\link[bacon]{bacon}} inflation
#' @param max_inflation (int) upper threshold values for \code{\link[bacon]{bacon}} inflation
#' @param decimal_places (int) lower threshold for average number of decimal places
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom stats quantile IQR
#'
#' @export
print_flags <- function(data_set, num_NA= 1000, min_beta=-2, max_beta=2, max_se=4, min_inflation=0.8, max_inflation=1.2, decimal_places=4,
                        save_path, gender, phenotype, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  #sets number of decimal places to 12 digits
  options(digits=12)

  #extract cohort names from data set
  cohort_names = as.character(unique(data_set$Cohort))
  full_path = paste0(save_path,phenotype,"_",gender,"_flags.txt")
  #terminal and log file output
  text = paste0("printing: flag file to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  sink(file= full_path, append= TRUE)
  #create flag file intro containing information on the tested parameters
  writeLines("Flag parameters tested:")
  writeLines("	- Different number of NAs in BETA, Pval and SE")
  writeLines(paste0("	- Number of NAs above ", num_NA))
  writeLines(paste0("	- min(BETA) below ", min_beta, " or max(BETA) above ", max_beta))
  writeLines(paste0("	- max(SE) above ", max_se))
  writeLines("	- more then 1 SE or BETA value above 75th + 1.5*IQR or below 25th percentile - 1.5*IQR")
  writeLines(paste0("	- Bacon inflation below ", min_inflation, " or above ", max_inflation))
  writeLines("	- number of cites of one chromosome = 0")
  writeLines("	- negative standard errors or p values")
  writeLines(c(paste0("	- average number of decimal places of BETA, Pval and SE smaller then ", decimal_places),"",""))

  for(i in 1:length(cohort_names))
  {
    #extract all the values from a cohort
    temp_data = data_set[data_set$Cohort == cohort_names[i],]

    ## test for NAs
    BETA_NA = sum(is.na(temp_data$BETA))
    PVAL_NA = sum(is.na(temp_data$P_VAL))
    SE_NA = sum(is.na(temp_data$SE))
    if((2*BETA_NA - PVAL_NA - SE_NA) != 0)
    {
      writeLines(paste0(cohort_names[i], " have different number of NA's. BETA: ", BETA_NA, " Pval: ", PVAL_NA, " SE: ", SE_NA))
    }
    if(BETA_NA > num_NA)
    {
      writeLines(paste0(cohort_names[i], " has high number of NAs (above ", num_NA, "): ", BETA_NA))
    }

    ## extreme BETA values
    if(min(temp_data$BETA, na.rm=TRUE) < min_beta || max(temp_data$BETA, na.rm=TRUE) > max_beta)
    {
      writeLines(paste0(cohort_names[i], " might have high/low BETA max/min. (smaller then ", min_beta, " or higher then ", max_beta,
                        "). max= ", max(temp_data$BETA, na.rm=TRUE), " min= ", min(temp_data$BETA, na.rm=TRUE)))
    }

    high_cutoff = stats::quantile(temp_data$BETA, probs=0.75, na.rm=TRUE) + (1.5*stats::IQR(temp_data$BETA, na.rm=TRUE))
    low_cutoff = stats::quantile(temp_data$BETA, probs=0.25, na.rm=TRUE) - (1.5*stats::IQR(temp_data$BETA, na.rm=TRUE))
    number_extreme = length(which(temp_data$BETA > high_cutoff | temp_data$BETA < low_cutoff))
    if(number_extreme > 1)
    {
      writeLines(paste0(cohort_names[i], " has ", number_extreme, " BETA values above 75th + 1.5*IQR or below 25th percentile - 1.5*IQR [",
                        low_cutoff, ",", high_cutoff,"]"))
    }

    ## extreme SE values
    if(max(temp_data$SE, na.rm=TRUE) > max_se)
    {
      writeLines(paste0(cohort_names[i], " might have high SE max. (higher then ", max_se, "). max= ", max(temp_data$SE, na.rm=TRUE)))
    }

    high_cutoff = stats::quantile(temp_data$SE, probs=0.75, na.rm=TRUE) + 1.5*stats::IQR(temp_data$SE, na.rm=TRUE)
    low_cutoff = stats::quantile(temp_data$SE, probs=0.25,  na.rm=TRUE) - 1.5*stats::IQR(temp_data$SE, na.rm=TRUE)
    number_extreme = length(which(temp_data$SE > high_cutoff | temp_data$SE < low_cutoff))
    if(number_extreme > 1)
    {
      writeLines(paste0(cohort_names[i], " has ", number_extreme, " SE values above 75th + 1.5*IQR or below 25th percentile - 1.5*IQR [",
                        low_cutoff, ",", high_cutoff,"]"))
    }

    ## inflation
    #inflation = qchisq(1-temp_data$P_VAL, 1)
    #inflation = median(inflation, na.rm= TRUE)/qchisq(0.5,1)
    #if (inflation > max_inflation || inflation < min_inflation)
    #{
      #writeLines(paste0(cohort_names[i], " has a high/low inflation value (above ", max_inflation, " or below ", min_inflation, "). inflation = ", inflation))
    #}

    se = temp_data$SE[!is.na(temp_data$BETA)]
    beta = temp_data$BETA[!is.na(temp_data$BETA)]
    beta = beta[!is.na(se)]
    se = se[!is.na(se)]

    inflation = bacon::inflation(bacon::bacon(NULL, beta, se))
    if(is.na(inflation))
    {
      writeLines(paste0(cohort_names[i], " BAChas produced NA's. Something is wrong."))
    }
    else
    {
      if(inflation > max_inflation || inflation < min_inflation)
      {
        writeLines(paste0(cohort_names[i], " has a high/low BACON inflation value (above ", max_inflation, " or below ", min_inflation, "). inflation = ", inflation))
      }
    }

    ## check for missing chromosomes
    chr_data = temp_data[!is.na(temp_data$BETA),]
    chr_data = chr_data[!is.na(chr_data$SE),]
    chr_data$CHR = as.numeric(chr_data$CHR)
    if(min(tabulate(chr_data$CHR, nbins=24)) == 0)
    {
      writeLines(paste0(cohort_names[i], " at least one chromosome is missing. Check chromosome plot."))
    }
    rm(chr_data)


    ## check decimal places
    if (mean(nchar(gsub("(.*\\.)|([0]*$)", "", as.character(temp_data$BETA))), na.rm=TRUE) < decimal_places)
    {
      writeLines(paste0(cohort_names[i], " average numbers of decimal places for BETA smaller then ", decimal_places, ". Average number of decimal places: ",
                        mean(nchar(gsub("(.*\\.)|([0]*$)", "", as.character(temp_data$BETA))), na.rm=TRUE)))
    }
    if (mean(nchar(gsub("(.*\\.)|([0]*$)", "", as.character(temp_data$P_VAL))), na.rm=TRUE) < decimal_places)
    {
      writeLines(paste0(cohort_names[i], " average numbers of decimal places for P_VAL smaller then ", decimal_places, ". Average number of decimal places: ",
                        mean(nchar(gsub("(.*\\.)|([0]*$)", "", as.character(temp_data$P_VAL))), na.rm=TRUE)))
    }
    if (mean(nchar(gsub("(.*\\.)|([0]*$)", "", as.character(temp_data$SE))), na.rm=TRUE) < decimal_places)
    {
      writeLines(paste0(cohort_names[i], " average numbers of decimal places for SE smaller then ", decimal_places, ". Average number of decimal places: ",
                        mean(nchar(gsub("(.*\\.)|([0]*$)", "", as.character(temp_data$SE))), na.rm=TRUE)))
    }

    ##check for negative p values
    if(all(temp_data$P_VAL > 0, na.rm=TRUE) == FALSE) { writeLines(paste0(cohort_names[i], " has negative p values")) }

    ##check for negative SE
    if(all(temp_data$SE > 0, na.rm=TRUE) == FALSE) { writeLines(paste0(cohort_names[i], " has negative SE values")) }

    writeLines(c("",""))
  }
  sink(file=NULL)
}
