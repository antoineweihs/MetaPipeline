#' @title Beta plots
#' @description  This function creates plots of the effect distributions. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param same_scale (bool) TRUE: the separate plots are all drawn using the same axis scales FALSE: the separate plots are all drawn using different scales
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param trim (bool) TRUE: will trim the density plots at local min and max FALSE: will plot density curve from global min to global max
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_line ggtitle ggsave facet_wrap aes
#'
#' @export
beta_plot <- function(data_set, same_scale=TRUE, save_path, gender, phenotype, trim=FALSE, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  num_cohorts = length(levels(data_set$Cohort))
  data_set$Cohort = as.character(data_set$Cohort)

  #all together
  text = "Plotting: Beta together"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  together = ggplot2::ggplot(data_set, ggplot2::aes(x=BETA, color=Cohort)) + ggplot2::geom_line(stat="density", trim=trim)
  together = together + ggplot2::ggtitle(paste0(phenotype, " ", gender, " BETA all together"))

  full_path = paste0(save_path,"BETA_together_",phenotype,"_",gender,".png")
  text = paste0("Saving: Beta together to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=together, width=594, height=420, units="mm", dpi="retina")

  #separate
  if(same_scale)
  {
    text = "Plotting: Beta separate (same scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=BETA)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", gender, " BETA separate"))

    full_path = paste0(save_path,"BETA_separate_samescale_",phenotype,"_",gender,".png")
    text = paste0("Saving: Beta separate (same scale) to: ", full_path)
  }
  if(!same_scale)
  {
    text = "Plotting: Beta separate (different scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=BETA)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., scales="free", ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", gender, " BETA separate DIFFERENT SCALES"))
    full_path = paste0(save_path,"BETA_separate_diffscale_",phenotype,"_",gender,".png")
    text = paste0("Saving: Beta separate (different scale) to: ", full_path)
  }

  ggplot2::ggsave(full_path, plot=separate, width=594, height=420, units="mm", dpi="retina")
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
}


#' @title Pvalue plots
#' @description  This function creates plots of the p-value distributions. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param same_scale (bool) True: the separate plots are all drawn using the same axis scales False: the separate plots are all drawn using different scales
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param trim (bool) TRUE: will trim the density plots at local min and max FALSE: will plot density curve from global min to global max
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#'
#' @importFrom ggplot2 ggplot geom_line ggtitle ggsave facet_wrap aes
#'
#' @export
pval_plot <- function(data_set,same_scale=TRUE, save_path, gender, phenotype, trim=FALSE, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  num_cohorts = length(levels(data_set$Cohort))
  data_set$Cohort = as.character(data_set$Cohort)

  #all together
  text = "Plotting: PVAL together"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  together = ggplot2::ggplot(data_set, ggplot2::aes(x=P_VAL, color=Cohort)) + ggplot2::geom_line(stat="density", trim=trim)
  together = together + ggplot2::ggtitle(paste0(phenotype, " ", gender, " PVAL all together"))
  full_path = paste0(save_path,"PVAL_together_",phenotype,"_",gender,".png")
  text = paste0("Saving: PVAL together to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=together, width=594, height=420, units="mm", dpi="retina")

  #separate
  if(same_scale)
  {
    text = "Plotting: PVAL separate (same scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=P_VAL)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., ncol = ceiling(sqrt(num_cohorts)))
    separate= separate + ggplot2::ggtitle(paste0(phenotype, " ", gender, " PVAL separate"))

    full_path = paste0(save_path,"PVAL_separate_samescale_",phenotype,"_",gender,".png")
    text = paste0("Saving: PVAL separate (same scale) to: ", full_path)
  }
  if(!same_scale)
  {
    text = "Plotting: PVAL separate (different scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=P_VAL)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., scales="free", ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", gender, " PVAL separate DIFFERENT SCALES"))

    full_path = paste0(save_path,"PVAL_separate_samescale_",phenotype,"_",gender,".png")
    text = paste0("Saving: PVAL separate (different scale) to: ", full_path)
  }

  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=separate, width=594, height=420, units="mm", dpi="retina")
}


#' @title Standard error plots
#' @description  This function creates plots of the standard error distributions. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param same_scale (bool) True: the separate plots are all drawn using the same axis scales False: the separate plots are all drawn using different scales
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param trim (bool) TRUE: will trim the density plots at local min and max FALSE: will plot density curve from global min to global max
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_line ggtitle ggsave facet_wrap aes
#'
#' @export
se_plot <- function(data_set,same_scale=TRUE, save_path, gender, phenotype, trim=FALSE, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  num_cohorts = length(levels(data_set$Cohort))
  data_set$Cohort = as.character(data_set$Cohort)

  #all together
  text = "Plotting: SE together"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  together = ggplot2::ggplot(data_set, ggplot2::aes(x=SE, color=Cohort)) + ggplot2::geom_line(stat="density", trim=trim)
  together = together + ggplot2::ggtitle(paste0(phenotype, " ", gender, " SE all together"))

  full_path = paste0(save_path,"SE_together_",phenotype,"_",gender,".png")
  text = paste0("Saving: SE together to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=together, width=594, height=420, units="mm", dpi="retina")

  #separate
  if(same_scale)
  {
    text = "Plotting: SE separate (same scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=SE)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", gender, " SE separate"))

    full_path = paste0(save_path,"SE_separate_samescale_",phenotype,"_",gender,".png")
    text = paste0("Saving: SE separate (same scale scale) to: ", full_path)
  }
  if(!same_scale)
  {
    text = "Plotting: SE separate (separate scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=SE)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., scales="free", ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", gender, " SE separate DIFFERENT SCALES"))

    full_path = paste0(save_path,"SE_separate_diffescale_",phenotype,"_",gender,".png")
    text = paste0("Saving: SE separate (separate scale) to: ", full_path)
  }

  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=separate, width=594, height=420, units="mm", dpi="retina")
}


#' @title Chromosome plots
#' @description This function creates plots visualising and quantifying the CpG site locations within the chromosomes. Used for the pre-meta quality
#'              control analysis.
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous ggtitle geom_text geom_bar ylab ggsave aes
#' @importFrom gridExtra grid.arrange
#'
#' @export
chromosome_plot <- function(data_set,save_path, gender, phenotype, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  fourfiftyk = data.frame(CHR=seq(1:24), Actual=c(46857,34810,25159,20464,24327,36611,30017,20950,9861,24388,28794,24539,12285,15078,15259,21969,27879,
                                                  5922,25521,10379,4243,8552,11232,416))
  EPIC = data.frame(CHR=seq(1:24), Actual=c(82013,64828,48896,36771,44720,54401,47560,38452,26167,42126,48894,44623,21040,29550,28741,37939,44435,14899,
                                            38550,22960,10300,18367,19090,537))
  cohort_names = as.character(unique(data_set$Cohort))
  for (i in 1:length(cohort_names))
  {
    text = paste0("Plotting: ", cohort_names[i], " chromosomes")
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

    temp_data = data_set[data_set$Cohort == cohort_names[i],]
    temp_data = temp_data[!is.na(temp_data$BETA),]
    temp_data = temp_data[!is.na(temp_data$SE),]
    temp_data$CHR = as.numeric(temp_data$CHR)
    temp_data$POS = as.numeric(temp_data$POS)
    if(all(temp_data$Array == "EPIC")) {isEPIC = TRUE}
    else if (all(temp_data$Array == "450")) {isEPIC = FALSE}
    else { stop(paste0("Unknown Array type", table(temp_data$Array))) }

    plot1 = ggplot2::ggplot(temp_data, ggplot2::aes(x=CHR, y=POS)) + ggplot2::geom_point() + ggplot2::scale_x_continuous(breaks=1:24, labels=1:24)
    plot1 = plot1 + ggplot2::ggtitle(paste0(phenotype, " ", gender, " ", cohort_names[i], " chromosomes"))
    if(isEPIC)
    {
      frequencies = data.frame(table(temp_data$CHR))
      names(frequencies) = c("CHR", "FREQ")
      frequencies = merge(EPIC, frequencies, by="CHR", all.x = TRUE)
      frequencies[is.na(frequencies)] = 0
      frequencies$label = paste0(frequencies$FREQ, "(",frequencies$FREQ - frequencies$Actual,")")
      plot2 = ggplot2::ggplot(frequencies) + ggplot2::geom_text(ggplot2::aes(x=CHR, y=Actual+1000, label=label))
      plot2 = plot2 + ggplot2::geom_bar(ggplot2::aes(x=CHR, y=Actual), stat="identity", color="red", fill="grey88")
      plot2 = plot2 + ggplot2::geom_bar(ggplot2::aes(x=CHR, y=FREQ), stat="identity", fill="black")
      plot2 = plot2 + ggplot2::scale_x_continuous(breaks=1:24, labels=1:24) + ggplot2::ylab("Frequency")
    }
    else
    {
      frequencies = data.frame(table(temp_data$CHR))
      names(frequencies) = c("CHR", "FREQ")
      frequencies = merge(fourfiftyk, frequencies, by="CHR", all.x = TRUE)
      frequencies[is.na(frequencies)] = 0
      frequencies$label = paste0(frequencies$FREQ, "(",frequencies$FREQ - frequencies$Actual,")")
      plot2 = ggplot2::ggplot(frequencies) + ggplot2::geom_text(ggplot2::aes(x=CHR, y=Actual+1000, label=label))
      plot2 = plot2 + ggplot2::geom_bar(ggplot2::aes(x=CHR, y=Actual), stat="identity", color="red", fill="grey88")
      plot2 = plot2 + ggplot2::geom_bar(ggplot2::aes(x=CHR, y=FREQ), stat="identity", fill="black")
      plot2 = plot2 + ggplot2::scale_x_continuous(breaks=1:24, labels=1:24) + ggplot2::ylab("Frequency")
    }

    plot = gridExtra::grid.arrange(plot1, plot2, nrow=2)

    full_path = paste0(save_path,"chromosomes_",cohort_names[i],"_",phenotype,"_",gender,".png")
    text = paste0("Saving: ", cohort_names[i], " chromosomes to: ", full_path)
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    ggplot2::ggsave(full_path, plot=plot, width=594, height=420, units="mm", dpi="retina")
  }
}


#' @title SE vs Size plots
#' @description This function creates plots the median cohort size versus the median cohort standard error. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <uni.antoine.weihs@@gmail.com>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param save_path (string) place where the outputs are saved
#' @param gender (string) gender currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text labs ggsave aes
#' @importFrom stats median
#'
#' @export
se_vs_size_plot <- function(data_set,save_path, gender, phenotype, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  cohort_names = as.character(unique(data_set$Cohort))
  median_se = rep(NA, length(cohort_names))
  median_size = rep(NA, length(cohort_names))
  for (i in 1:length(cohort_names))
  {
    median_size[i] = stats::median(data_set[data_set$Cohort == cohort_names[i],]$Size, na.rm=TRUE)
    median_se[i] = stats::median(data_set[data_set$Cohort == cohort_names[i],]$SE, na.rm=TRUE)
  }

  text = paste0("Plotting: 1 / standard error vs size")
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  temp_data = data.frame(median_SE = median_se, median_Size= 1/median_size, label=cohort_names)
  plot = ggplot2::ggplot(temp_data, ggplot2::aes(x=median_SE, y=median_Size, label=label)) + ggplot2::geom_point()
  plot = plot + ggplot2::ggtitle(paste0(phenotype, " ", gender, " ", " se vs size")) + ggplot2::geom_text(hjust=0, nudge_x = 0.05)
  plot = plot + ggplot2::labs(x="median standard error", y="1 / median size")

  full_path = paste0(save_path,"se_vs_size_",phenotype,"_",gender,".png")
  text = paste0("Saving plot to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=plot, width=594, height=420, units="mm", dpi="retina")
}
