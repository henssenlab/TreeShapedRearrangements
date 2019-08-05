rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(circlize)
library(BSgenome.Hsapiens.UCSC.hg19)

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Circles.Rdata")

plot_onesample_onesvcaller_plus = function(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, fname){
  
  hg19sub = keepStandardChromosomes(GRangesForBSGenome("hg19"), pruning.mode = "coarse")
  bins <- tileGenome(seqinfo(hg19sub), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
  
  mycols = rand_color(n=100, luminosity = "dark")
  bgcolor = "gray98"
  
  palmtrees_of_interest = palmtrees %>% filter(Sample == sample, SVCaller == svcaller)
  
  pdf(file=fname)
  circos.clear()
  circos.par("start.degree" = 90)
  
  # Introduce gap between chrY and chr1 at the top to label tracks manually in Illustrator
  gapsize = 10 # intr
  circos.par("start.degree" = 90-gapsize/2, "gap.after" = c(rep(1,23), gapsize))
  
  circos.initializeWithIdeogram(species = "hg19")
  #circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
  text(0, 0, "", cex = 1) 
  
  # Plotting Copy Number Variants
  ascat_cnv$logFC = log10((ascat_cnv$TumorTotalCopyNumber + 0.001) / (ascat_cnv$NormalTotalCopyNumber + 0.001))
  ascat_cnv_bed = ascat_cnv %>% filter(Sample == sample) %>% dplyr::select(Chr, Start, End, logFC)
  ascat_cnv_bed_gr = makeGRangesFromDataFrame(ascat_cnv_bed, keep.extra.columns = T, seqinfo = seqinfo(hg19sub))
  ascat_cnv_bed_gr = binnedAverage(bins=bins, coverage(ascat_cnv_bed_gr, weight=ascat_cnv_bed_gr$logFC), "y") %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  
  circos.genomicTrack(as.data.frame(ascat_cnv_bed_gr), 
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value[[1]],  ytop.column = 1, ybottom = 0,
                                           col = ifelse(value[[1]] > 0, "firebrick4", ifelse(value[[1]] < 0, "steelblue4", "lightgray")),
                                           border = ifelse(value[[1]] > 0, "firebrick4", ifelse(value[[1]] < 0, "steelblue4", "lightgray")), ...)
                      },
                      track.height=0.5*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  #circos.text(sector.index="chr1",track.index = 1,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
  #            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "CNV", facing = "clockwise", 
  #            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
  
  # Plotting WGS eccDNA
  wgs_eccDNA_bed = circles %>% filter(Sample == sample, Method  == "WGSCircles", CircleClass == "eccDNA") %>%
    dplyr::select(Chr, Start, End) %>% as.data.frame() %>%
    mutate(Chr = gsub("chrMT", "chrM", Chr))
  wgs_eccDNA_bed_gr = wgs_eccDNA_bed %>% makeGRangesFromDataFrame(., seqinfo=seqinfo(hg19sub))
  wgs_eccDNA_bed_gr_binned = binnedAverage(bins=bins, coverage(wgs_eccDNA_bed_gr), "y")  %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  circos.genomicTrack(wgs_eccDNA_bed_gr_binned,
                      panel.fun = function(region, value, ...){
                        circos.genomicLines(region, value, col="darkolivegreen4", ...)
                      },
                      track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  
  # Plotting WGS ecDNA
  wgs_ecDNA_bed = circles %>% filter(Sample == sample, Method == "WGSCircles", CircleClass == "ecDNA") %>%
    dplyr::select(Chr, Start, End) %>% as.data.frame()
  circos.genomicTrack(wgs_ecDNA_bed,
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col="darkorange3", border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  
  if (circles %>% filter(Sample == sample, Method  == "CircleSeq") %>% nrow() > 0){
  # Plotting Circle-seq eccDNA
  circleseq_eccDNA_bed = circles %>% filter(Sample == sample, Method  == "CircleSeq", CircleClass == "eccDNA") %>%
    dplyr::select(Chr, Start, End) %>% as.data.frame()
  circleseq_eccDNA_bed_gr = circleseq_eccDNA_bed %>% makeGRangesFromDataFrame(., seqinfo=seqinfo(hg19sub))
  circleseq_eccDNA_bed_gr_binned = binnedAverage(bins=bins, coverage(circleseq_eccDNA_bed_gr), "y")  %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  circos.genomicTrack(circleseq_eccDNA_bed_gr_binned,
                      panel.fun = function(region, value, ...){
                        circos.genomicLines(region, value, col="darkorchid4", ...)
                      },
                      track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  # Plotting Circle-seq ecDNA
  circleseq_ecDNA_bed = circles %>% filter(Sample == sample, Method == "CircleSeq", CircleClass == "ecDNA") %>%
    dplyr::select(Chr, Start, End) %>% as.data.frame()
  circos.genomicTrack(circleseq_ecDNA_bed,
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col="deeppink3", border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  }
  
  # Plot Palm Trees
  palmtrees_bed = palmtrees_of_interest[,c("Chr", "FirstElement", "LastElement")]
  
  # this is just for plotting
  diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
  if (sum(diff<5000000)>0){
    palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff)/2
    palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff)/2
  }
  
  if (nrow(palmtrees_bed)>0){
    palmtrees_bed$value1 = 1:nrow(palmtrees_of_interest)
    palmtrees_bed$value2 = palmtrees_of_interest$PalmTreeID
  }else{
    palmtrees_bed$value1 = data.frame()
    palmtrees_bed$value2 = data.frame()
  }
  colnames(palmtrees_bed) = c("chr", "start", "end", "value1", "value2")
  
  #bed_palmtree_A = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
  #bed_palmtree_B = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
  
  circos.genomicTrack(as.data.frame(palmtrees_bed),
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col="indianred4", border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  
  this_tx_circles = tx_circles %>% filter(Sample == sample, SVCaller == svcaller, 
                                          CircleMethod == circlemethod)
  
  bed_tx_A = this_tx_circles %>% mutate(PosACopy = PosA) %>% dplyr::select(ChrA, PosA, PosACopy, CircleGenomeClass)
  bed_tx_A = bed_tx_A %>% mutate(CircleGenomeClass = 
                                   ifelse(CircleGenomeClass == "genome-genome", "steelblue4",
                                          ifelse(CircleGenomeClass == "circle-genome", "firebrick4", 
                                                 ifelse(CircleGenomeClass == "circle-circle", "green3", "lightgray"))))
  colnames(bed_tx_A) = c("chr", "start", "end", "value1")
  bed_tx_B = this_tx_circles %>% mutate(PosBCopy = PosB) %>% dplyr::select(ChrB, PosB, PosBCopy, CircleGenomeClass)
  colnames(bed_tx_B) = c("chr", "start", "end", "value1")
  circos.genomicLink(as.data.frame(bed_tx_A), as.data.frame(bed_tx_B), col=bed_tx_A$value1)
  
  title(paste(sample, svcaller))
  dev.off()
  circos.clear()
}

sample = "NB2013"
svcaller = "Union"
circlemethod = "CircleSeq" # or WGSCircles
txdouble_ptinfo = txptinfo
plot_onesample_onesvcaller_plus(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPro/", svcaller, "/", circlemethod, "/", svcaller, "_", circlemethod, "_", sample, "_Circos.pdf"))

sample = "NBL27"
svcaller = "Union"
circlemethod = "WGSCircles" 
txdouble_ptinfo = txptinfo
plot_onesample_onesvcaller_plus(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPro/", svcaller, "/", circlemethod, "/", svcaller, "_", circlemethod, "_", sample, "_Circos.pdf"))

sample = "NBL04"
svcaller = "Union"
circlemethod = "WGSCircles" 
txdouble_ptinfo = txptinfo
plot_onesample_onesvcaller_plus(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPro/", svcaller, "/", circlemethod, "/", svcaller, "_", circlemethod, "_", sample, "_Circos.pdf"))


# 
# samples = unique(tx_circles$Sample)
# for (sample in samples){
#   svcallers = tx_circles %>% filter(Sample == sample) %>% .$SVCaller %>% unique()
#   for (svcaller in svcallers){
#     circlemethods = tx_circles %>% filter(Sample == sample, SVCaller == svcaller) %>% .$CircleMethod %>% unique()
#     for (circlemethod in circlemethods){
#       plot_onesample_onesvcaller_plus(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPro/", svcaller, "/", circlemethod, "/", svcaller, "_", circlemethod, "_", sample, "_Circos.pdf"))
#     }
#   }
# }


# samples = unique(as.character(tx$Sample))
# svcallers = unique(as.character(tx$SVCaller))
# for (sai in 1:length(samples)){
#   print(samples[sai])
#   for (svi in 1:length(svcallers)){
#     if (nrow(tx %>% filter(Sample == samples[sai], SVCaller == svcallers[svi]))>0){
#       plot_onesample_onesvcaller_plus(circles, ascat_cnv, palmtrees, txptinfo, txdouble, samples[sai], svcallers[svi], paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPlus/", as.character(svcallers[svi]), "/", as.character(samples[sai]), as.character(svcallers[svi]), "PalmTreeCircosPlus.pdf"))
#     }
#   }
# }

circles %>% filter(Sample == "NB2013", Method=="WGSCircles", CircleClass =="ecDNA") %>% View()

# for debugging:
# circles %>% filter(Sample == "NB2013", Method=="CircleSeq", CircleClass =="eccDNA") %>% 
#   ggplot(aes(x=Start, y=abs(Start-End)))+
#   geom_point() + 
#   facet_wrap(Chr ~.)
####


plot_onesample_onesvcaller_plus_color_tx_by_palmtree = function(circles, ascat_cnv, palmtrees, txdouble_ptinfo, txdouble, sample, svcaller, fname){
  
  hg19sub = keepStandardChromosomes(GRangesForBSGenome("hg19"), pruning.mode = "coarse")
  bins <- tileGenome(seqinfo(hg19sub), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
  
  mycols = rand_color(n=100, luminosity = "dark")
  
  palmtrees_of_interest = palmtrees %>% filter(Sample == sample, SVCaller == svcaller)
  
  pdf(file=fname)
  circos.clear()
  circos.par("start.degree" = 90)
  
  # Introduce gap between chrY and chr1 at the top to label tracks manually in Illustrator
  gapsize = 30 # intr
  circos.par("start.degree" = 90-gapsize/2, "gap.after" = c(rep(1,23), gapsize))
  
  #circos.initializeWithIdeogram(species = "hg19")
  circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
  text(0, 0, "", cex = 1) 
  
  # Plotting Copy Number Variants
  ascat_cnv$logFC = log10((ascat_cnv$TumorTotalCopyNumber + 0.001) / (ascat_cnv$NormalTotalCopyNumber + 0.001))
  ascat_cnv_bed = ascat_cnv %>% filter(Sample == sample) %>% dplyr::select(Chr, Start, End, logFC)
  ascat_cnv_bed_gr = makeGRangesFromDataFrame(ascat_cnv_bed, keep.extra.columns = T, seqinfo = seqinfo(hg19sub))
  ascat_cnv_bed_gr = binnedAverage(bins=bins, coverage(ascat_cnv_bed_gr, weight=ascat_cnv_bed_gr$logFC), "y") %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  
  circos.genomicTrack(as.data.frame(ascat_cnv_bed_gr), 
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value[[1]],  ytop.column = 1, ybottom = 0,
                                           col = ifelse(value[[1]] > 0, "#e41a1c", ifelse(value[[1]] < 0, "#377eb8", "lightgray")),
                                           border = ifelse(value[[1]] > 0, "#e41a1c", ifelse(value[[1]] < 0, "#377eb8", "lightgray")), ...)
                      },
                      track.height=0.5*circos.par("track.height"))
  
  #circos.text(sector.index="chr1",track.index = 1,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
  #            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "CNV", facing = "clockwise", 
  #            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
  
  # # Plotting WGS eccDNA
  # wgs_eccDNA_bed = circles %>% filter(Sample == sample, Method  == "WGSCircles", CircleClass == "eccDNA") %>%
  #   dplyr::select(Chr, Start, End) %>% as.data.frame() %>%
  #   mutate(Chr = gsub("chrMT", "chrM", Chr))
  # wgs_eccDNA_bed_gr = wgs_eccDNA_bed %>% makeGRangesFromDataFrame(., seqinfo=seqinfo(hg19sub))
  # wgs_eccDNA_bed_gr_binned = binnedAverage(bins=bins, coverage(wgs_eccDNA_bed_gr), "y")  %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  # circos.genomicTrack(wgs_eccDNA_bed_gr_binned,
  #                     panel.fun = function(region, value, ...){
  #                       circos.genomicLines(region, value, col="#4daf4a", ...)
  #                     },
  #                     track.height=0.25*circos.par("track.height"))
  
  # Plotting WGS ecDNA
  wgs_ecDNA_bed = circles %>% filter(Sample == sample, Method == "WGSCircles", CircleClass == "ecDNA") %>%
    dplyr::select(Chr, Start, End) %>% as.data.frame()
  circos.genomicTrack(wgs_ecDNA_bed,
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col="#4daf4a", border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"))
  
  # Plotting Circle-seq eccDNA
  if (circles %>% filter(Sample == sample, Method  == "CircleSeq", CircleClass == "eccDNA") %>% nrow() > 0){
    # circleseq_eccDNA_bed = circles %>% filter(Sample == sample, Method  == "CircleSeq", CircleClass == "eccDNA") %>%
    #   dplyr::select(Chr, Start, End) %>% as.data.frame()
    # circleseq_eccDNA_bed_gr = circleseq_eccDNA_bed %>% makeGRangesFromDataFrame(., seqinfo=seqinfo(hg19sub))
    # circleseq_eccDNA_bed_gr_binned = binnedAverage(bins=bins, coverage(circleseq_eccDNA_bed_gr), "y")  %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
    # circos.genomicTrack(circleseq_eccDNA_bed_gr_binned,
    #                     panel.fun = function(region, value, ...){
    #                       circos.genomicLines(region, value, col="#984ea3", ...)
    #                     },
    #                     track.height=0.25*circos.par("track.height"))
    
    # Plotting Circle-seq ecDNA
    circleseq_ecDNA_bed = circles %>% filter(Sample == sample, Method == "CircleSeq", CircleClass == "ecDNA") %>%
      dplyr::select(Chr, Start, End) %>% as.data.frame()
    circos.genomicTrack(circleseq_ecDNA_bed,
                        panel.fun = function(region, value, ...){
                          circos.genomicRect(region, value, col="#984ea3", border=NA, ...)
                        },
                        ylim=c(0,0.2), track.height=0.25*circos.par("track.height"))
  }
  
  # Plot Palm Trees
  palmtrees_bed = palmtrees_of_interest[,c("Chr", "FirstElement", "LastElement")]
  
  # this is just for plotting
  diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
  if (sum(diff<5000000)>0){
    palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff)/2
    palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff)/2
  }
  
  if (nrow(palmtrees_bed)>0){
    palmtrees_bed$value1 = 1:nrow(palmtrees_of_interest)
    palmtrees_bed$value2 = palmtrees_of_interest$PalmTreeID
  }else{
    palmtrees_bed$value1 = data.frame()
    palmtrees_bed$value2 = data.frame()
  }
  colnames(palmtrees_bed) = c("chr", "start", "end", "value1", "value2")
  
  #bed_palmtree_A = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
  #bed_palmtree_B = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
  
  bed_palmtree_A = data.frame()
  bed_palmtree_B = data.frame()
  if (nrow(palmtrees_bed) > 0){
    for (i in 1:nrow(palmtrees_bed)){
      palmtree_tx = txdouble_ptinfo %>% filter(PalmTreeID == as.character(palmtrees_bed[i, "value2"]))
      palmtree_tx$PalmTreeIndex = i
      palmtree_tx = palmtree_tx %>% dplyr::select(PalmTreeChrom, PalmTreePos, TargetChrom, TargetPos, PalmTreeIndex) %>% distinct()
      temp = palmtree_tx[,c("PalmTreeChrom", "PalmTreePos", "PalmTreePos", "PalmTreeIndex")]
      colnames(temp) = c("chr", "start", "end", "value1")
      bed_palmtree_A = rbind(bed_palmtree_A, temp)
      colnames(bed_palmtree_A) = c("chr", "start", "end", "value1")
      temp = palmtree_tx[,c("TargetChrom", "TargetPos", "TargetPos", "PalmTreeIndex")]
      colnames(temp) = c("chr", "start", "end", "value1")
      bed_palmtree_B = rbind(bed_palmtree_B, temp)
      colnames(bed_palmtree_B) = c("chr", "start", "end", "value1")
    }
  }
  
  circos.genomicTrack(as.data.frame(palmtrees_bed),
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col=mycols[value[[1]]], border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"))
  
  nopalmtree_tx = tx_original %>% filter(Sample == sample, SVCaller == svcaller, isdefchrom(ChrA), isdefchrom(ChrB))
  bed_nopalmtree_A = nopalmtree_tx[,c("ChrA", "PosA", "PosA")]
  if (nrow(bed_nopalmtree_A) > 0){
    bed_nopalmtree_A$value1 = 0
    colnames(bed_nopalmtree_A) = c("chr", "start", "end", "value1")
  }
  
  bed_nopalmtree_B = nopalmtree_tx[,c("ChrB", "PosB", "PosB")]
  if (nrow(bed_nopalmtree_B) > 0){
    bed_nopalmtree_B$value1 = 0
    colnames(bed_nopalmtree_B) = c("chr", "start", "end", "value1")
  }
  circos.genomicLink(bed_nopalmtree_A, bed_nopalmtree_B, col="lightgray")
  
  if (nrow(palmtrees_bed) > 0) circos.genomicLink(bed_palmtree_A, bed_palmtree_B, col=mycols[bed_palmtree_A[,"value1"]])
  
  title(paste(sample, svcaller))
  dev.off()
  circos.clear()

  }

###############


plot_onesample_onesvcaller_plus2 = function(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, fname){
  
  hg19sub = keepStandardChromosomes(GRangesForBSGenome("hg19"), pruning.mode = "coarse")
  bins <- tileGenome(seqinfo(hg19sub), tilewidth=1000000, cut.last.tile.in.chrom=TRUE)
  
  mycols = rand_color(n=100, luminosity = "dark")
  bgcolor = "gray98"
  
  palmtrees_of_interest = palmtrees %>% filter(Sample == sample, SVCaller == svcaller)
  
  pdf(file=fname)
  circos.clear()
  circos.par("start.degree" = 90)
  
  # Introduce gap between chrY and chr1 at the top to label tracks manually in Illustrator
  gapsize = 10 # intr
  circos.par("start.degree" = 90-gapsize/2, "gap.after" = c(rep(1,23), gapsize))
  
  #circos.initializeWithIdeogram(species = "hg19")
  circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
  text(0, 0, "", cex = 1) 
  
  # Plotting Copy Number Variants
  ascat_cnv$logFC = log10((ascat_cnv$TumorTotalCopyNumber + 0.001) / (ascat_cnv$NormalTotalCopyNumber + 0.001))
  ascat_cnv_bed = ascat_cnv %>% filter(Sample == sample) %>% dplyr::select(Chr, Start, End, logFC)
  ascat_cnv_bed_gr = makeGRangesFromDataFrame(ascat_cnv_bed, keep.extra.columns = T, seqinfo = seqinfo(hg19sub))
  ascat_cnv_bed_gr = binnedAverage(bins=bins, coverage(ascat_cnv_bed_gr, weight=ascat_cnv_bed_gr$logFC), "y") %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  
  circos.genomicTrack(as.data.frame(ascat_cnv_bed_gr), 
                      panel.fun = function(region, value, ...) {
                        circos.genomicRect(region, value[[1]],  ytop.column = 1, ybottom = 0,
                                           col = ifelse(value[[1]] > 0, "firebrick", ifelse(value[[1]] < 0, "steelblue", "lightgray")),
                                           border = ifelse(value[[1]] > 0, "firebrick", ifelse(value[[1]] < 0, "steelblue", "lightgray")), ...)
                      },
                      track.height=0.5*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  #circos.text(sector.index="chr1",track.index = 1,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
  #            get.cell.meta.data("cell.ylim")-max(get.cell.meta.data("cell.ylim"))/2, labels = "CNV", facing = "clockwise", 
  #            niceFacing = TRUE, adj = c(0,0),cex = 0.5)
  
  # Plotting WGS eccDNA
  # wgs_eccDNA_bed = circles %>% filter(Sample == sample, Method  == "WGSCircles", CircleClass == "eccDNA") %>%
  #   dplyr::select(Chr, Start, End) %>% as.data.frame() %>%
  #   mutate(Chr = gsub("chrMT", "chrM", Chr))
  # wgs_eccDNA_bed_gr = wgs_eccDNA_bed %>% makeGRangesFromDataFrame(., seqinfo=seqinfo(hg19sub))
  # wgs_eccDNA_bed_gr_binned = binnedAverage(bins=bins, coverage(wgs_eccDNA_bed_gr), "y")  %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
  # circos.genomicTrack(wgs_eccDNA_bed_gr_binned,
  #                     panel.fun = function(region, value, ...){
  #                       circos.genomicLines(region, value, col="darkolivegreen4", ...)
  #                     },
  #                     track.height=0.25*circos.par("track.height"),
  #                     bg.border=NA, bg.col = bgcolor)
  
  
  # Plotting WGS ecDNA
  wgs_ecDNA_bed = circles %>% filter(Sample == sample, Method == "WGSCircles", CircleClass == "ecDNA") %>%
    dplyr::select(Chr, Start, End) %>% as.data.frame()
  circos.genomicTrack(wgs_ecDNA_bed,
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col="darkorange3", border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  
  if (circles %>% filter(Sample == sample, Method  == "CircleSeq") %>% nrow() > 0){
    # Plotting Circle-seq eccDNA
    # circleseq_eccDNA_bed = circles %>% filter(Sample == sample, Method  == "CircleSeq", CircleClass == "eccDNA") %>%
    #   dplyr::select(Chr, Start, End) %>% as.data.frame()
    # circleseq_eccDNA_bed_gr = circleseq_eccDNA_bed %>% makeGRangesFromDataFrame(., seqinfo=seqinfo(hg19sub))
    # circleseq_eccDNA_bed_gr_binned = binnedAverage(bins=bins, coverage(circleseq_eccDNA_bed_gr), "y")  %>% as.data.frame() %>% dplyr::select(seqnames, start, end, y)
    # circos.genomicTrack(circleseq_eccDNA_bed_gr_binned,
    #                     panel.fun = function(region, value, ...){
    #                       circos.genomicLines(region, value, col="darkorchid4", ...)
    #                     },
    #                     track.height=0.25*circos.par("track.height"),
    #                     bg.border=NA, bg.col = bgcolor)
    
    # Plotting Circle-seq ecDNA
    circleseq_ecDNA_bed = circles %>% filter(Sample == sample, Method == "CircleSeq", CircleClass == "ecDNA") %>%
      dplyr::select(Chr, Start, End) %>% as.data.frame()
    circos.genomicTrack(circleseq_ecDNA_bed,
                        panel.fun = function(region, value, ...){
                          circos.genomicRect(region, value, col="deeppink3", border=NA, ...)
                        },
                        ylim=c(0,0.2), track.height=0.25*circos.par("track.height"),
                        bg.border=NA, bg.col = bgcolor)
  }
  
  # Plot Palm Trees
  palmtrees_bed = palmtrees_of_interest[,c("Chr", "FirstElement", "LastElement")]
  
  # this is just for plotting
  diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
  if (sum(diff<5000000)>0){
    palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff)/2
    palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff)/2
  }
  
  if (nrow(palmtrees_bed)>0){
    palmtrees_bed$value1 = 1:nrow(palmtrees_of_interest)
    palmtrees_bed$value2 = palmtrees_of_interest$PalmTreeID
  }else{
    palmtrees_bed$value1 = data.frame()
    palmtrees_bed$value2 = data.frame()
  }
  colnames(palmtrees_bed) = c("chr", "start", "end", "value1", "value2")
  
  #bed_palmtree_A = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
  #bed_palmtree_B = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
  
  circos.genomicTrack(as.data.frame(palmtrees_bed),
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col="indianred4", border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.25*circos.par("track.height"),
                      bg.border=NA, bg.col = bgcolor)
  
  
  this_tx_circles = tx_circles %>% filter(Sample == sample, SVCaller == svcaller, 
                                          CircleMethod == circlemethod)
  
  bed_tx_A = this_tx_circles %>% mutate(PosACopy = PosA) %>% dplyr::select(ChrA, PosA, PosACopy, CircleGenomeClass)
  bed_tx_A = bed_tx_A %>% mutate(CircleGenomeClass = 
                                   ifelse(CircleGenomeClass == "genome-genome", "steelblue",
                                          ifelse(CircleGenomeClass == "circle-genome", "firebrick", 
                                                 ifelse(CircleGenomeClass == "circle-circle", "green2", "lightgray"))))
  colnames(bed_tx_A) = c("chr", "start", "end", "value1")
  bed_tx_B = this_tx_circles %>% mutate(PosBCopy = PosB) %>% dplyr::select(ChrB, PosB, PosBCopy, CircleGenomeClass)
  colnames(bed_tx_B) = c("chr", "start", "end", "value1")
  circos.genomicLink(as.data.frame(bed_tx_A), as.data.frame(bed_tx_B), col=bed_tx_A$value1)
  
  title(paste(sample, svcaller))
  dev.off()
  circos.clear()
}

sample = "NB2013"
svcaller = "Union"
circlemethod = "CircleSeq" # or WGSCircles
txdouble_ptinfo = txptinfo
plot_onesample_onesvcaller_plus2(circles, ascat_cnv, palmtrees, tx_circles, sample, svcaller, circlemethod, paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPro/", svcaller, "/", circlemethod, "/", svcaller, "_", circlemethod, "_", sample, "_Circos_Cleaner.pdf"))

