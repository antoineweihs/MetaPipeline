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
  num_cohorts = length(unique(data_set$Cohort))
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
  save(together, file = paste0(full_path, ".RData"))
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

  save(separate, file = paste0(full_path, ".RData"))
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
  num_cohorts = length(unique(data_set$Cohort))
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
  
  save(together, file = paste0(full_path, ".RData"))
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

  save(separate, file = paste0(full_path, ".RData"))
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
  num_cohorts = length(unique(data_set$Cohort))
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

  save(together, file = paste0(full_path, ".RData"))
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

  save(separate, file = paste0(full_path, ".RData"))
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
    temp_data$MAPINFO = as.numeric(temp_data$MAPINFO)
    if(all(temp_data$Array == "EPIC")) {isEPIC = TRUE}
    else if (all(temp_data$Array == "450")) {isEPIC = FALSE}
    else { stop(paste0("Unknown Array type", table(temp_data$Array))) }

    plot1 = ggplot2::ggplot(temp_data, ggplot2::aes(x=CHR, y=MAPINFO)) + ggplot2::geom_point() + ggplot2::scale_x_continuous(breaks=1:24, labels=1:24)
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

    save(plot, file = paste0(full_path, ".RData"))
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
    median_size[i] = sqrt(stats::median(data_set[data_set$Cohort == cohort_names[i],]$Size, na.rm=TRUE)) # use sqrt(N)
    median_se[i] = stats::median(data_set[data_set$Cohort == cohort_names[i],]$SE, na.rm=TRUE)
  }

  text = paste0("Plotting: 1 / standard error vs size")
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  temp_data = data.frame(median_SE = 1/median_se, median_Size= median_size, label=cohort_names)
  plot = ggplot2::ggplot(temp_data, ggplot2::aes(x=median_SE, y=median_Size, label=label)) + ggplot2::geom_point()
  plot = plot + ggplot2::ggtitle(paste0(phenotype, " ", stratum, " ", " weight vs size")) + ggplot2::geom_text(hjust=0, vjust=0, check_overlap = TRUE)
  plot = plot + ggplot2::labs(x="1/median standard error", y="sqrt(median size)")

  full_path = paste0(save_path,"se_vs_size_",phenotype,"_",stratum,".png")
  text = paste0("Saving plot to: ", full_path)
  if(verbose) {writeLines(text)}
  if(print_log) {cat(text, file=log_path, append=TRUE, sep="\n")}

  save(plot, file = paste0(full_path, ".RData"))
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
double_manhattan <- function(x, chr="CHR", bp="MAPINFO", p="P", markername="MARKERNAME", beta="BETA",
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
  d = d[!is.na(d$CHR),]
  d = d[!is.na(d$BP),]

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

  if (save_plot) 
  {
    save(myplot, file = paste0(save_path,"_manhattan_plot.png.RData"))
    ggplot2::ggsave(paste0(save_path,"_manhattan_plot.png"), plot=myplot, width=15.5, height=8.61, units="cm", dpi="retina")
  }

  return(0)
}


#' @title annotation plots for CpG sites
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @description plotting function that creates annotation plots of CpG sites using the UCSC database
#'
#' @param result (data.frame) data set created by \code{\link{meta}}
#' @param id (string) CpG ID of the CpG site that should be plotted
#' @param phenotype (string) name of the phenotype analysed
#' @param width (int) width of the window size around the CpG site that should be plotted
#' @param FDR (bool) TRUE: uses FDR correction FALSE: won't
#' @param verbose (bool) TRUE: will print output to terminal FALSE: won't
#' @param print_log (bool) TRUE: prints logs to log_path FALSE: won't
#' @param log_path (string) full path + file name of log file
#' @param save_dest (string) path to where the plots should be saved
#'
#' @importFrom Gviz IdeogramTrack AnnotationTrack UcscTrack displayPars HighlightTrack DataTrack plotTracks
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom rtracklayer browserSession genome ucscTableQuery GRangesForUCSCGenome getTable
#' @importFrom IRanges IRanges ranges
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors mcols
#' @export
annotation_plot <- function(result, id, phenotype, width=50000, FDR =F,
                            verbose=T, print_log = F, log_path = "log.txt",
                            save_dest = "C:/Users/weihsa/Documents/")
{
  "%nin%" = Negate("%in%")
  if("MAPINFO" %nin% names(result))
  {
    result = merge(result, Masterfile, by="Markername")
  }

  if(FDR & ("FDR" %nin% names(result)))
  {
    result = result[!is.na(result$Pval_phenotype),]
    result$FDR = stats::p.adjust(result$Pval_phenotype, "BH")
    text = paste0("FDR missing from data set. Calculating FDR values")
    if(verbose) writeLines(text)
    if(print_log) cat(text, file=log_path, append=TRUE, sep="\n")
  }

  text = paste0("Plotting: ", phenotype, " ", id)
  if(verbose) writeLines(text)
  if(print_log) cat(text, file=log_path, append=TRUE, sep="\n")

  #get location of cpg site of interest
  chr = result$CHR[result$Markername == id]
  location = result$MAPINFO[result$Markername == id]

  #extract cpg sites from result data set which lie +- 50,000bp from the site of interest
  temp = result[result$CHR == chr,]
  temp = temp[(temp$MAPINFO >= (location - width) & temp$MAPINFO <= (location + width)),]


  #reduce data set in order to fit the expected format
  if(FDR)
  {
    temp = temp[, c(which(names(temp) == "Markername"), which(names(temp) == "CHR"),
                    which(names(temp) == "MAPINFO"), which(names(temp) == "FDR"))]
  }
  else
  {
    temp = temp[, c(which(names(temp) == "Markername"), which(names(temp) == "CHR"),
                  which(names(temp) == "MAPINFO"), which(names(temp) == "Pval_phenotype"))]
  }

  temp = temp[order(temp$MAPINFO),]
  names(temp) = c("TargetID", "CHR", "MAPINFO", "Pval")

  gen = "hg19"
  coord = as.numeric(as.character(temp$MAPINFO))
  #little chromosome plot at the top
  itrack <- Gviz::IdeogramTrack(genome = gen, chromosome = chr)


  #pval plot
  pval = -log10(as.numeric(as.character(temp$Pval)))
  mygroups = rep(1, length(temp$Pval))
  mygroups[temp$TargetID == id] = 2
  if(FDR)
  {
    cpgTrack = Gviz::DataTrack(data = pval, start = coord, width=0, chromosome = chr, genome = gen, name = "-log10(FDR)", group=mygroups,
                               baseline = -log10(1.1e-07), col.baseline = "#FF2D00", col = "black", type="p", symbol=temp$TargetID,
                               cex=1.2, pch=21, bg="slateblue3")
  }
  else
  {
    cpgTrack = Gviz::DataTrack(data = pval, start = coord, width=0, chromosome = chr, genome = gen, name = "-log10(P-Values)", group=mygroups,
                               baseline = -log10(1.1e-07), col.baseline = "#FF2D00", col = "black", type="p", symbol=temp$TargetID,
                               cex=1.2, pch=21, bg="slateblue3")
  }

  #ucsc gene plot
  mystart = (min(coord)-width)
  myend = (max(coord)+width)
  knownGenes = Gviz::UcscTrack(genome=gen, chromosome=chr, table ="ncbiRefSeq", track = 'NCBI RefSeq', from=mystart, to=myend,
                               trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                               symbol="name", transcript="name", strand="strand", fill="#8282d2", name="RefSeq Genes",
                               showID=T, geneSymbol=T, stacking = 'pack')
  z = IRanges::ranges(knownGenes)
  S4Vectors::mcols(z)$symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, gsub("\\.[1-9]$", "", S4Vectors::mcols(z)$symbol), "SYMBOL","REFSEQ")
  IRanges::ranges(knownGenes) = z

  #open session
  mySession = rtracklayer::browserSession()
  rtracklayer::genome(mySession) = gen
  #cpg Islands
  cpgIslandquery = rtracklayer::ucscTableQuery(mySession, "CpG Islands", rtracklayer::GRangesForUCSCGenome(gen, paste0("chr",chr), IRanges::IRanges(mystart, myend)))
  cpgIslandTable = rtracklayer::getTable(cpgIslandquery)
  if(dim(cpgIslandTable)[1] == 0) cpgIsland = Gviz::AnnotationTrack(NULL, name = "CpG Islands")
  else
  {
    cpgIslandGranges = GenomicRanges::makeGRangesFromDataFrame(cpgIslandTable, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
    cpgIsland = Gviz::AnnotationTrack(cpgIslandGranges, name = "CpG Islands")
  }

  #chromHMM
  #chromHMMquery = rtracklayer::ucscTableQuery(mySession, "Broad ChromHMM", rtracklayer::GRangesForUCSCGenome(gen, paste0("chr", chr), IRanges::IRanges(mystart, myend)))
  #chromHMMTable = rtracklayer::getTable(chromHMMquery)
  #chromHMMTable$hex = sapply(strsplit(as.character(chromHMMTable$itemRgb), ","), function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  #chromHMMGranges = GenomicRanges::makeGRangesFromDataFrame(chromHMMTable, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
  #chromHMM = Gviz::AnnotationTrack(chromHMMGranges, name="  Broad ChromHMM", fill=chromHMMTable$hex, stacking = "dense", col.line = "black", collapse=F)

  #Common SNPs (build 151)
  SNPquery = rtracklayer::ucscTableQuery(mySession,  "Common SNPs(151)", rtracklayer::GRangesForUCSCGenome(gen, paste0("chr", chr), IRanges::IRanges(mystart, myend)))
  SNPtable = rtracklayer::getTable(SNPquery)
  SNPGranges = GenomicRanges::makeGRangesFromDataFrame(SNPtable, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
  SNP = Gviz::AnnotationTrack(SNPGranges, name="Common SNPs (151)", stacking = "dense", col = NULL, fill = "black")


  Gviz::displayPars(cpgIsland) = list(cex.title = 0.5, rotation.title = 0)
  #Gviz::displayPars(chromHMM) = list(cex.title = 0.5, rotation.title = 0)
  Gviz::displayPars(SNP) = list(cex.title = 0.5, rotation.title = 0)
  Gviz::displayPars(cpgTrack) = list(cex.title = 0.5)
  Gviz::displayPars(knownGenes) = list(cex.title = 0.5, rotation.title = 0)
  #ht <- Gviz::HighlightTrack(trackList = list(cpgTrack, knownGenes, cpgIsland, chromHMM, SNP), start = location, end = location, chromosome = chr, col = "black")
  ht <- Gviz::HighlightTrack(trackList = list(cpgTrack, knownGenes, cpgIsland, SNP), start = location, end = location, chromosome = chr, col = "black")

  #actual plot
  tracklist = list(itrack, ht)
  from = min(coord)-width
  to = max(coord)+width
  transAnnot = "symbol"
  save(tracklist, from, to, transAnnot, file = paste0(save_dest, phenotype, "_", id, ".png.RData"))
  png(file = paste0(save_dest, phenotype, "_", id, ".png"), res=320, width=16, height=16, units = "cm")
  Gviz::plotTracks(tracklist, from = from, to = to, transcriptAnnotation=transAnnot)
  dev.off()
}
