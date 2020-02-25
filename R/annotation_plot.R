#' @title annotation plots for CpG sites
#' @author Antoine Weihs <antoine.weihs@@uni-greifswald.de>
#' @description plotting function that creates annotation plots of CpG sites using the UCSC database
#'
#' @param result (data.frame) data set created by \code{\link{meta}}
#' @param id (string) CpG ID of the CpG site that should be plotted
#' @param phenotype (string) name of the phenotype analysed
#' @param width (int) width of the window size around the CpG site that should be plotted
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
annotation_plot <- function(result, id, phenotype, width=50000,
                            verbose=T, print_log = F, log_path = "log.txt",
                            save_dest = "C:/Users/weihsa/Documents/")
{
  "%nin%" = Negate("%in%")
  if("MAPINFO" %nin% names(result)) result = merge(result, Masterfile, by="Markername")

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
  temp = temp[, c(which(names(temp) == "Markername"), which(names(temp) == "CHR"),
                  which(names(temp) == "MAPINFO"), which(names(temp) == "Pval_phenotype"))]
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
  cpgTrack = Gviz::DataTrack(data = pval, start = coord, width=0, chromosome = chr, genome = gen, name = "-log10(P-Values)", group=mygroups,
                             baseline = -log10(1.1e-07), col.baseline = "#FF2D00", col = c("blue", "red"), type="p", symbol=temp$TargetID)

  #ucsc gene plot
  mystart = (min(coord)-width)
  myend = (max(coord)+width)
  knownGenes = Gviz::UcscTrack(genome=gen, chromosome=chr, table ="ncbiRefSeq", track = 'NCBI RefSeq', from=mystart, to=myend,
                               trackType="GeneRegionTrack", rstarts="exonStarts", rends="exonEnds", gene="name",
                               symbol="name", transcript="name", strand="strand", fill="#8282d2", name="UCSC Genes",
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
  cpgIslandGranges = GenomicRanges::makeGRangesFromDataFrame(cpgIslandTable, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
  cpgIsland = Gviz::AnnotationTrack(cpgIslandGranges, name = "CpG Islands")

  #chromHMM
  chromHMMquery = rtracklayer::ucscTableQuery(mySession, "Broad ChromHMM", rtracklayer::GRangesForUCSCGenome(gen, paste0("chr", chr), IRanges::IRanges(mystart, myend)))
  chromHMMTable = rtracklayer::getTable(chromHMMquery)
  chromHMMTable$hex = sapply(strsplit(as.character(chromHMMTable$itemRgb), ","), function(x) rgb(x[1], x[2], x[3], maxColorValue = 255))
  chromHMMGranges = GenomicRanges::makeGRangesFromDataFrame(chromHMMTable, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
  chromHMM = Gviz::AnnotationTrack(chromHMMGranges, name="  Broad ChromHMM", fill=chromHMMTable$hex, stacking = "dense", col.line = "black", collapse=F)

  #SNPs (build 153)
  SNPquery = rtracklayer::ucscTableQuery(mySession, "dbSNP 153", rtracklayer::GRangesForUCSCGenome(gen, paste0("chr", chr), IRanges::IRanges(mystart, myend)))
  SNPtable = rtracklayer::getTable(SNPquery)
  SNPGranges = GenomicRanges::makeGRangesFromDataFrame(SNPtable, keep.extra.columns = T, start.field = "chromStart", end.field = "chromEnd")
  SNP = Gviz::AnnotationTrack(SNPGranges, name="  dbSNP 153", stacking = "dense", col = NULL, fill = "black")


  Gviz::displayPars(cpgIsland) = list(cex.title = 0.5, rotation.title = 0)
  Gviz::displayPars(chromHMM) = list(cex.title = 0.5, rotation.title = 0)
  Gviz::displayPars(SNP) = list(cex.title = 0.5, rotation.title = 0)
  Gviz::displayPars(cpgTrack) = list(cex.title = 0.5)
  Gviz::displayPars(knownGenes) = list(cex.title = 0.5, rotation.title = 0)
  ht <- Gviz::HighlightTrack(trackList = list(cpgTrack, knownGenes, cpgIsland, chromHMM, SNP), start = location, end = location, chromosome = chr, col = "black")

  #actual plot
  png(file = paste0(save_dest, phenotype, "_", id, ".png"), res=320, width=16, height=16, units = "cm")
  Gviz::plotTracks(list(itrack, ht), from = min(coord)-width, to = max(coord)+width, transcriptAnnotation="symbol")
  dev.off()
}
