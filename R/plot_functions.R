#' @title Beta plots
#' @description  This function creates plots of the effect distributions. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param same_scale (bool) TRUE: the separate plots are all drawn using the same axis scales FALSE: the separate plots are all drawn using different scales
#' @param save_path (string) place where the outputs are saved
#' @param stratum (string) stratum currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param trim (bool) TRUE: will trim the density plots at local min and max FALSE: will plot density curve from global min to global max
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_line ggtitle ggsave facet_wrap aes
#'
#' @export
beta_plot <- function(data_set, same_scale=TRUE, save_path, stratum, phenotype, trim=FALSE, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  num_cohorts = length(levels(data_set$Cohort))
  data_set$Cohort = as.character(data_set$Cohort)

  #all together
  text = "Plotting: Beta together"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  together = ggplot2::ggplot(data_set, ggplot2::aes(x=BETA, color=Cohort)) + ggplot2::geom_line(stat="density", trim=trim)
  together = together + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " BETA all together"))

  full_path = paste0(save_path,"BETA_together_",phenotype,"_",stratum,".png")
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
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " BETA separate"))

    full_path = paste0(save_path,"BETA_separate_samescale_",phenotype,"_",stratum,".png")
    text = paste0("Saving: Beta separate (same scale) to: ", full_path)
  }
  if(!same_scale)
  {
    text = "Plotting: Beta separate (different scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=BETA)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., scales="free", ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " BETA separate DIFFERENT SCALES"))
    full_path = paste0(save_path,"BETA_separate_diffscale_",phenotype,"_",stratum,".png")
    text = paste0("Saving: Beta separate (different scale) to: ", full_path)
  }

  ggplot2::ggsave(full_path, plot=separate, width=594, height=420, units="mm", dpi="retina")
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
}


#' @title Pvalue plots
#' @description  This function creates plots of the p-value distributions. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param same_scale (bool) True: the separate plots are all drawn using the same axis scales False: the separate plots are all drawn using different scales
#' @param save_path (string) place where the outputs are saved
#' @param stratum (string) stratum currently analysed
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
pval_plot <- function(data_set,same_scale=TRUE, save_path, stratum, phenotype, trim=FALSE, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  num_cohorts = length(levels(data_set$Cohort))
  data_set$Cohort = as.character(data_set$Cohort)

  #all together
  text = "Plotting: PVAL together"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  together = ggplot2::ggplot(data_set, ggplot2::aes(x=P_VAL, color=Cohort)) + ggplot2::geom_line(stat="density", trim=trim)
  together = together + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " PVAL all together"))
  full_path = paste0(save_path,"PVAL_together_",phenotype,"_",stratum,".png")
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
    separate= separate + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " PVAL separate"))

    full_path = paste0(save_path,"PVAL_separate_samescale_",phenotype,"_",stratum,".png")
    text = paste0("Saving: PVAL separate (same scale) to: ", full_path)
  }
  if(!same_scale)
  {
    text = "Plotting: PVAL separate (different scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=P_VAL)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., scales="free", ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " PVAL separate DIFFERENT SCALES"))

    full_path = paste0(save_path,"PVAL_separate_samescale_",phenotype,"_",stratum,".png")
    text = paste0("Saving: PVAL separate (different scale) to: ", full_path)
  }

  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=separate, width=594, height=420, units="mm", dpi="retina")
}


#' @title Standard error plots
#' @description  This function creates plots of the standard error distributions. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param same_scale (bool) True: the separate plots are all drawn using the same axis scales False: the separate plots are all drawn using different scales
#' @param save_path (string) place where the outputs are saved
#' @param stratum (string) stratum currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param trim (bool) TRUE: will trim the density plots at local min and max FALSE: will plot density curve from global min to global max
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_line ggtitle ggsave facet_wrap aes
#'
#' @export
se_plot <- function(data_set,same_scale=TRUE, save_path, stratum, phenotype, trim=FALSE, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
{
  num_cohorts = length(levels(data_set$Cohort))
  data_set$Cohort = as.character(data_set$Cohort)

  #all together
  text = "Plotting: SE together"
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  together = ggplot2::ggplot(data_set, ggplot2::aes(x=SE, color=Cohort)) + ggplot2::geom_line(stat="density", trim=trim)
  together = together + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " SE all together"))

  full_path = paste0(save_path,"SE_together_",phenotype,"_",stratum,".png")
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
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " SE separate"))

    full_path = paste0(save_path,"SE_separate_samescale_",phenotype,"_",stratum,".png")
    text = paste0("Saving: SE separate (same scale scale) to: ", full_path)
  }
  if(!same_scale)
  {
    text = "Plotting: SE separate (separate scale)"
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    separate = ggplot2::ggplot(data_set, ggplot2::aes(x=SE)) + ggplot2::geom_line(stat="density", trim=trim)
    separate = separate + ggplot2::facet_wrap(Cohort ~ ., scales="free", ncol = ceiling(sqrt(num_cohorts)))
    separate = separate + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " SE separate DIFFERENT SCALES"))

    full_path = paste0(save_path,"SE_separate_diffescale_",phenotype,"_",stratum,".png")
    text = paste0("Saving: SE separate (separate scale) to: ", full_path)
  }

  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=separate, width=594, height=420, units="mm", dpi="retina")
}


#' @title Chromosome plots
#' @description This function creates plots visualising and quantifying the CpG site locations within the chromosomes. Used for the pre-meta quality
#'              control analysis.
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param save_path (string) place where the outputs are saved
#' @param stratum (string) stratum currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_point scale_x_continuous ggtitle geom_text geom_bar ylab ggsave aes
#' @importFrom gridExtra grid.arrange
#'
#' @export
chromosome_plot <- function(data_set,save_path, stratum, phenotype, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
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
    plot1 = plot1 + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " ", cohort_names[i], " chromosomes"))
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

    full_path = paste0(save_path,"chromosomes_",cohort_names[i],"_",phenotype,"_",stratum,".png")
    text = paste0("Saving: ", cohort_names[i], " chromosomes to: ", full_path)
    if(verbose) {writeLines(text)}
    if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
    ggplot2::ggsave(full_path, plot=plot, width=594, height=420, units="mm", dpi="retina")
  }
}


#' @title SE vs Size plots
#' @description This function creates plots the median cohort size versus the median cohort standard error. Used for the pre-meta quality control analysis.
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @param data_set (data.frame) data set created by \code{\link{load_files}}
#' @param save_path (string) place where the outputs are saved
#' @param stratum (string) stratum currently analysed
#' @param phenotype (string) phenotype currently analysed
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log
#'
#' @importFrom ggplot2 ggplot geom_point ggtitle geom_text labs ggsave aes
#' @importFrom stats median
#'
#' @export
se_vs_size_plot <- function(data_set,save_path, stratum, phenotype, verbose=TRUE, print_log=FALSE, log_path="./log.txt")
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
  plot = plot + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " ", " se vs size")) + ggplot2::geom_text(hjust=0, nudge_x = 0.05)
  plot = plot + ggplot2::labs(x="median standard error", y="1 / median size")

  full_path = paste0(save_path,"se_vs_size_",phenotype,"_",stratum,".png")
  text = paste0("Saving plot to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}
  ggplot2::ggsave(full_path, plot=plot, width=594, height=420, units="mm", dpi="retina")
}


#' @title double Manhattan Plot
#' @description This function creates a double manahattan plot for \code{\link{post_process}}
#'              The base code was extracted from the manhattan() function from the qqman r package
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @references Stephen D. Turner, 2018, "qqman: an R package for visualizing GWAS results usinf Q-Q and manhattan plots", Jorunal of Open Source Software, 3(25), 731, doi:10:21105/joss.00731
#' @param x (data.frame) input data set
#' @param chr (string) name of the coloumn in \code{x} that contains the chromosome numbers
#' @param bp (string) name of the coloumn in \code{x} that contains the genome position
#' @param p (string) name of the coloumn in \code{x} that contains the p-values
#' @param markername (string) name of the column in \code{x} that contains the cpg site names
#' @param col (string vector) name of the two colours the manhattan plot shall be coloured in
#' @param FDRcorr (bool) TRUE: will perform an FDR correction on the p-values and use these FALSE: will use the p-values
#' @param cutoff (double) significance level
#' @param strict_cutoff (double) strict significance level
#' @param title (string) title of the manhattan plot
#' @param logp (bool) TRUE will use a logarithmic scale on the y-axis FALSE: won't
#' @param marktop (bool) TRUE: will label the cpg sites above the cutoff FALSE: won't
#' @param save_plot (bool) TRUE: saves the plot to \code{save_path} FALSE: won't
#' @param save_path (string) save path for the manhattan plot
#'
#' @importFrom ggplot2 ggplot geom_point scale_color_manual geom_text xlab ggsave aes ggtitle scale_x_continuous geom_hline theme element_blank ylab
#'
#' @export
double_manhattan <- function(x, chr="CHR", bp="BP", p="P", markername="MARKERNAME", beta="BETA",
                             col=c("gray10", "gray60"), FDRcorr = TRUE,
                             cutoff=0.05, strict_cutoff=0.001, title = "Manhattan plot",
                             logp=TRUE, marktop = FALSE, save_plot=TRUE, save_path="./",...)
{

  temp <- function(x) return(as.numeric(as.character(x)))
  d = data.frame(CHR=temp(x[[chr]]), BP=temp(x[[bp]]), P=temp(x[[p]]), beta=temp(x[[beta]]), pos = rep(NA, length(x[[bp]])),
                 index = rep(NA, length(x[[bp]])), MARKERNAME=x[[markername]],  stringsAsFactors = FALSE)

  if(length(d$CHR) == 0)
  {
    return(1)
  }

  d = d[!is.na(d$beta),]
  d = d[!is.na(d$P),]

  if(FDRcorr) d$P = p.adjust(d$P, "BH")

  # Set positions, ticks, and labels for plotting
  ## Sort and keep only values where is numeric.
  d <- d[order(d$CHR, d$BP), ]
  if (logp) {
    d$P = d$P + 1E-323  ##kann ich das machen? (falls p wert = 0)
    d$logp <- -log10(d$P)
    if (!is.null(cutoff)) cutoff = -log10(cutoff)
    if (!is.null(strict_cutoff)) strict_cutoff = -log10(strict_cutoff)
  } else {
    d$logp <- d$P
  }

  d$logp = d$logp * sign(d$beta)

  ind = 0
  for(i in unique(d$CHR))
  {
    ind= ind+1
    d[d$CHR == i,]$index = ind
  }

  # This section sets up positions and ticks. Ticks should be placed in the
  # middle of a chromosome. The a new pos column is added that keeps a running
  # sum of the positions of each successive chromsome. For example:
  # chr bp pos
  # 1   1  1
  # 1   2  2
  # 2   1  3
  # 2   2  4
  # 3   1  5
  nchr = length(unique(d$CHR))
  if (nchr==1) { ## For a single chromosome
    d$pos=d$BP
    xlabel = paste('Chromosome',unique(d$CHR),'position')
  } else { ## For multiple chromosomes
    lastbase=0
    ticks=NULL
    for (i in unique(d$index)) {
      if (i==1) {
        d[d$index==i, ]$pos=d[d$index==i, ]$BP
      } else {
        ## chromosome position maybe not start at 1, eg. 9999. So gaps may be produced.
        lastbase = lastbase +max(d[d$index==(i-1),"BP"])
        d[d$index == i,"BP"] = d[d$index == i,"BP"]-min(d[d$index==i,"BP"]) +1
        d[d$index == i, "pos"] = d[d$index == i,"BP"] + lastbase

      }
    }
    ticks <-tapply(d$pos,d$index,quantile,probs=0.5)
    xlabel = 'Chromosome'
    labs <- unique(d$CHR)
  }

  ## create plot
  if(marktop & !is.null(cutoff))
  {
    myplot = ggplot2::ggplot(d, ggplot2::aes(x=pos, y=logp, label=MARKERNAME)) +  ggplot2::geom_point(ggplot2::aes(color=as.factor(CHR)), shape=20, size=0.25)
    myplot = myplot + ggplot2::scale_color_manual(values=rep(col, max(labs))) + ggplot2::xlab("Chromosome")
    myplot = myplot + ggplot2::geom_text(ggplot2::aes(label=ifelse(abs(logp) > cutoff, as.character(MARKERNAME), "")), hjust=0, vjust=0)
  }
  if(!marktop)
  {
    myplot = ggplot2::ggplot(d, ggplot2::aes(x=pos, y=logp)) +  ggplot2::geom_point(ggplot2::aes(color=as.factor(CHR)), shape=".")
    myplot = myplot + ggplot2::scale_color_manual(values=rep(col, max(labs))) + ggplot2::xlab("Chromosome")

  }

  myplot = myplot + ggplot2::scale_x_continuous(label=labs, breaks=ticks) + ggplot2::geom_hline(yintercept=0, color="white", alpha=0.2)
  myplot = myplot + ggplot2::theme(legend.position="none", panel.border=ggplot2::element_blank(), panel.grid.major.x=ggplot2::element_blank(),
                                   panel.grid.minor.x=ggplot2::element_blank())
  myplot = myplot + ggplot2::ggtitle(title)

  if (logp & FDRcorr) myplot = myplot + ggplot2::ylab("-log10(FDR) * effect direction")
  if (logp & !FDRcorr) myplot = myplot + ggplot2::ylab("-log10(p)  * effect direction")
  if (!logp & FDRcorr) myplot = myplot + ggplot2::ylab("FDR  * effect direction")
  if (!logp & !FDRcorr) myplot = myplot + ggplot2::ylab("p value  * effect direction")


  if (!is.null(cutoff)) myplot = myplot + ggplot2::geom_hline(yintercept=cutoff, col="blue") + ggplot2::geom_hline(yintercept=-cutoff, col="blue")
  if (!is.null(strict_cutoff)) myplot = myplot + ggplot2::geom_hline(yintercept=strict_cutoff, col="red") + ggplot2::geom_hline(yintercept=-strict_cutoff, col="red")

  if (save_plot) ggplot2::ggsave(paste0(save_path,"_manhattan_plot.png"), plot=myplot, width=15.5, height=8.61, units="cm", dpi="retina")

  return(0)
}
