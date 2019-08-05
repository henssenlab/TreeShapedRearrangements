#load("~/Desktop/PalmTrees/Results/Figures/RelationToWGSCircles/WGSCircleOverlap_RandomSampling500.Rdata")
#load("~/Desktop/PalmTrees/Results/Figures/RelationToCSCircles/CSCircleOverlap_RandomSampling500.Rdata")
#rm(list=setdiff(ls(), c("palmtrees_onlymatched", "txptinfo_onlymatched", "cscircles", "merged_wgscircles")))
rm(list=ls())
library(dplyr)
library(ggbio)
library(GenomicRanges)
library(biomaRt)
library(Homo.sapiens) # is based on hg19
data("ideoCyto", package = "biovizBase")
data(genesymbol, package = "biovizBase")
library(BSgenome.Hsapiens.UCSC.hg19)
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Circles.Rdata")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/MergedTx.Rdata")

circlegenomeclass = tx_circles %>%
  filter(CircleMethod == "CircleSeq") %>%
  dplyr::select(BPID, CircleGenomeClass)

plot_palm_tree = function(thispalmtree, thissample, thissvcaller, thischr, thisstart, thisend, palmtrees, tx_original, cscircles, merged_wgscircles, ascat_cnv, circlegenomeclass){
  
  # thispalmtree = palmtrees[[i, "PalmTreeID"]]
  # thissample = palmtrees[[i,"Sample"]]
  # thissvcaller = palmtrees[[i,"SVCaller"]]
  # thischr = palmtrees[[i,"Chr"]]
  # palmtreelength = palmtrees[[i,"LastElement"]] - palmtrees[[i,"FirstElement"]]
  # thisstart = palmtrees[[i,"FirstElement"]] - ceiling(0.25 * palmtreelength)
  # thisend = palmtrees[[i,"LastElement"]] +  ceiling(0.25 * palmtreelength)

  thispalmtree = "NB2013_Union_chr2:14550246-16220971"
  thissample = "NB2013"
  thissvcaller = "Union"
  thischr = "chr2"
  palmtreelength = 16220971 - 14550246
  thisstart = 14550246 - ceiling(0.25 * palmtreelength)
  thisend = 16220971 +  ceiling(0.25 * palmtreelength)

  # thispalmtree = "NB2050_Union_chr2:15740877-16158228"
  # thissample = "NB2050"
  # thissvcaller = "Union"
  # thischr = "chr2"
  # palmtreelength = 16158228 - 15740877
  # thisstart = 15740877 - ceiling(0.25 * palmtreelength)
  # thisend = 16158228 +  ceiling(0.25 * palmtreelength)

  #library(dplyr)
  #library(ggbio)
  #library(GenomicRanges)
  #library(Homo.sapiens) # is based on hg19
  #data("ideoCyto", package = "biovizBase")
  #data(genesymbol, package = "biovizBase")
  #library(BSgenome.Hsapiens.UCSC.hg19)
  #source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
  
  # i=3
  # thispalmtree=palmtrees[[i, "PalmTreeID"]]
  # thissample=palmtrees[[i,"Sample"]]
  # thissvcaller=palmtrees[[i,"SVCaller"]] 
  # thischr=palmtrees[[i,"Chr"]] 
  # thisstart=palmtrees[[i,"FirstElement"]] - ceiling(0.25 * palmtreelength) 
  # thisend=palmtrees[[i,"LastElement"]] +  ceiling(0.25 * palmtreelength) 

  
  thisstart = max(0, thisstart)
  thisend = min(thisend, seqlengths(BSgenome.Hsapiens.UCSC.hg19)[[thischr]])

  ### BUILD TRANSLOCATION TRACK
  
  tx_original_2 = tx_original
  colnames(tx_original_2)[4:7] = c("ChrB", "PosB", "ChrA", "PosA") 
  tx_original_double = rbind(tx_original, tx_original_2)
  tx_original_double = inner_join(tx_original_double, circlegenomeclass)
  tx_granges = tx_original_double %>% 
    filter(Sample==thissample, SVCaller ==thissvcaller, ChrA == thischr, hasOverlap(PosA, PosA, thisstart, thisend)) %>% 
    dplyr::select(ChrA, PosA, PosA, CircleGenomeClass) %>% distinct()
  if (nrow(tx_granges) > 0){
    tx_granges = makeGRangesFromDataFrame(tx_granges, seqnames.field = "ChrA", start.field="PosA", end.field="PosA", ignore.strand = T, keep.extra.columns = T)
    tx.track = autoplot(tx_granges, geom="rect", size=0.1, aes(fill=CircleGenomeClass, color=CircleGenomeClass), size=2) + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank()) + scale_fill_manual(values=c("circle-circle"="green4", "circle-genome"="firebrick2", "genome-genome"="steelblue")) + scale_color_manual(values=c("circle-circle"="green4", "circle-genome"="firebrick2", "genome-genome"="steelblue")) + guides(color=F, fill=F)
  } else {
    tx_granges = GRanges()
    tx.track = ggplot() + theme_null()
  }
  
  
  ### BUILD PALM TREE TRACK
  
  pt_granges = palmtrees %>% 
    filter(Sample==thissample, SVCaller==thissvcaller) %>% 
    filter(Sample==thissample, Chr == thischr, hasOverlap(FirstElement, LastElement, thisstart, thisend)) %>% 
    dplyr::select(Chr, FirstElement, LastElement)
  if (nrow(pt_granges)>0){
    pt_granges = makeGRangesFromDataFrame(pt_granges, seqnames.field = "Chr", start.field="FirstElement", end.field="LastElement", ignore.strand = T)
    pt.track = autoplot(pt_granges, geom="rect", color=NA, fill="gold") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
  } else {
    pt_granges = GRanges()
    pt.track = ggplot() + theme_null()    
  }

  ### BUILD CIRCLE-SEQ TRACK eccDNA
  circ_ecc_granges = cscircles_ecc %>%
    filter(Sample==thissample, CircleChr == thischr, hasOverlap(CircleStart, CircleEnd, thisstart, thisend)) %>% 
    dplyr::select(CircleChr, CircleStart, CircleEnd)
  if (nrow(circ_ecc_granges)>0){
    circ_ecc_granges = makeGRangesFromDataFrame(circ_ecc_granges, seqnames.field = "CircleChr", start.field="CircleStart", end.field="CircleEnd", ignore.strand = T)
    circ.ecc.track = autoplot(circ_ecc_granges, geom="rect", color=NA, fill="#984ea3") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
  } else {
    circ_ecc_granges = GRanges()
    circ.ecc.track = ggplot() + theme_null()
  }
  
  ### BUILD CIRCLE-SEQ TRACK ecDNA
  
  circ_ec_granges = cscircles_ec %>%
    filter(Sample==thissample, CircleChr == thischr, hasOverlap(CircleStart, CircleEnd, thisstart, thisend)) %>% 
    dplyr::select(CircleChr, CircleStart, CircleEnd)
  if (nrow(circ_ec_granges)>0){
    circ_ec_granges = makeGRangesFromDataFrame(circ_ec_granges, seqnames.field = "CircleChr", start.field="CircleStart", end.field="CircleEnd", ignore.strand = T)
    circ.ec.track = autoplot(circ_ec_granges, geom="rect", color=NA, fill="#984ea3") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
  } else {
    circ.ec.granges = GRanges()
    circ.ec.track = ggplot() + theme_null()
  }
  
  
  # seqlevels(circ_granges, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  # seqinfo(circ_granges) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  # circ_granges <- trim(circ_granges)
  
  
  ### BUILD WGS eccDNA TRACK
  wgscirc_ecc_granges = merged_wgscircles_ecc %>%
    filter(Sample==thissample, CircleChr == thischr, hasOverlap(CircleStart, CircleEnd, thisstart, thisend)) %>% 
    dplyr::select(CircleChr, CircleStart, CircleEnd)
  if (nrow(wgscirc_ecc_granges) > 0){
    wgscirc_ecc_granges = makeGRangesFromDataFrame(wgscirc_ecc_granges, seqnames.field = "CircleChr", start.field="CircleStart", end.field="CircleEnd", ignore.strand = T)
    wgscirc.ecc.track = autoplot(wgscirc_ecc_granges, geom="rect", color=NA, fill= "#4daf4a") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
  }else{
    wgscirc_ecc_granges = GRanges()
    wgscirc.ecc.track = ggplot() + theme_null()
  }
  
  ### BUILD WGS ecDNA TRACK
  wgscirc_ec_granges = merged_wgscircles_ec %>%
    filter(Sample==thissample, CircleChr == thischr, hasOverlap(CircleStart, CircleEnd, thisstart, thisend)) %>% 
    dplyr::select(CircleChr, CircleStart, CircleEnd)
  if (nrow(wgscirc_ec_granges) > 0){
    wgscirc_ec_granges = makeGRangesFromDataFrame(wgscirc_ec_granges, seqnames.field = "CircleChr", start.field="CircleStart", end.field="CircleEnd", ignore.strand = T)
    wgscirc.ec.track = autoplot(wgscirc_ec_granges, geom="rect", color=NA, fill= "#4daf4a") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
  }else{
    wgscirc_ec_granges = GRanges()
    wgscirc.ec.track = ggplot() + theme_null()
  }

  ### BUILD COPY NUMBER TRACK
  
  ascat_granges = ascat_cnv %>%
    filter(Sample==thissample, Chr == thischr, hasOverlap(Start, End, thisstart, thisend)) %>% 
    mutate(y = log10((TumorTotalCopyNumber + 0.01) / (NormalTotalCopyNumber + 0.01))) %>%
    dplyr::select(Chr, Start, End, y)
  if (nrow(ascat_granges)>0){
    ascat_granges = makeGRangesFromDataFrame(ascat_granges, seqnames.field = "Chr", start.field="Start", end.field="End", ignore.strand = T, keep.extra.columns=T)
    ascat.track = autoplot(ascat_granges, geom="bar", aes(y=y), fill="#e41a1c", color="#e41a1c", ylab="", size=0.1) + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
  } else {
    ascat_granges = GRanges()
    ascat.track = ggplot() + theme_null()
  }
  
  
  ### BUILD TRANSCRIPTS TRACK
  
  if (((thisend-thisstart)<=10000000) & ((thisend-thisstart)>20)){
    gene.track = try(autoplot(Homo.sapiens, 
                              which = GRanges(seqnames=thischr,
                                              ranges=IRanges(start = thisstart, end = thisend),
                                              strand="*"),
                              names.expr="SYMBOL",
                              label.color="black") +
                       xlim(thisstart,thisend) + 
                       theme_clear(),    
                     silent=TRUE)
    if (is.character(gene.track)){ # if the try caught an error 
      gene.track = ggplot(data=data.frame(x = c(thisstart, thisend), y=c(0,0)), aes(x=x, y=y)) + geom_point(alpha=0) + xlim(thisstart, thisend) + annotation_custom(grid::textGrob("There is probably no \n transcripts in the region of interest"), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_clear() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) 
    }
  } else {
    gene.track=ggplot(data=data.frame(x = c(thisstart, thisend), y=c(0,0)), aes(x=x, y=y)) + geom_point(alpha=0) + xlim(thisstart, thisend) + annotation_custom(grid::textGrob("Too small or \n too large region \n to retrieve transcripts."), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_clear() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
  }
  
  
  ### PUT ALL TRACKS TOGETHER
  
  if (thissample %in% cscircles_ecc$Sample){
    tracks("Translocations" = tx.track, "Palm Tree \n Region" = pt.track, "eccDNA \n (Circle-Seq)" = circ.ecc.track, "ecDNA \n (Circle-Seq)" = circ.ec.track, "eccDNA \n (WGS)" = wgscirc.ecc.track, "ecDNA \n (WGS)" = wgscirc.ec.track, "Copy Number \n logFC" = ascat.track, "Transcripts" = gene.track, 
           heights = c(1,1,1,1,1,1,1,2),
           track.bg.color=NULL,
           xlim=c(thisstart,thisend),
           xlab=thischr,
           title = paste0("Palm Tree Region ", thispalmtree),
           label.text.cex = 1,
           label.bg.fill = "white",
           label.bg.color = "white") +
      theme(axis.line=element_line(color="black")) +
      ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/StackFigures/", thispalmtree, "_StackFigure.pdf"), height=10, width=10)
  } else {
    tracks("Translocations" = tx.track, "Palm Tree \n Region" = pt.track, "eccDNA \n (WGS)" = wgscirc.ecc.track, "ecDNA \n (WGS)" = wgscirc.ec.track, "Copy Number \n logFC" = ascat.track, "Transcripts" = gene.track, 
           heights = c(1,1,1,1,1,2),
           track.bg.color=NULL,
           xlab=thischr,
           title = paste0("Palm Tree Region ", thispalmtree),
           label.text.cex = 1,
           label.bg.fill = "white",
           label.bg.color = "white") +
      theme(axis.line=element_line(color="black")) +
      ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/StackFigures/", thispalmtree, "_StackFigure.pdf"), height=10, width=10)
  }
  
}


## PLOT PALM TREE REGIONS
for (i in 1:nrow(palmtrees)){
  palmtreelength = palmtrees[[i,"LastElement"]] - palmtrees[[i,"FirstElement"]]
  plot_palm_tree(palmtrees[[i, "PalmTreeID"]],
                 palmtrees[[i,"Sample"]], 
                 palmtrees[[i,"SVCaller"]], 
                 palmtrees[[i,"Chr"]], 
                 palmtrees[[i,"FirstElement"]] - ceiling(0.25 * palmtreelength), 
                 palmtrees[[i,"LastElement"]] +  ceiling(0.25 * palmtreelength), 
                 palmtrees, tx_original, cscircles, merged_wgscircles, ascat_cnv, circlegenomeclass)
}


###################################


# plot_target_site = function(txptinfo, i, thispalmtree, thissample, thissvcaller, thischr, thisstart, thisend, palmtrees, tx_original, cscircles, merged_wgscircles, ascat_cnv){
#   
#   thispalmtree = txptinfo[i, "PalmTreeID"]
#   thissample = txptinfo[i, "Sample"]
#   thissvcaller = txptinfo[i, "SVCaller"]
#   thischr = txptinfo[i, "TargetChrom"]
#   thisstart = txptinfo[i, "TargetPos"]-500000
#   thisend = txptinfo[i, "TargetPos"]+500000
#   
#   
#   library(dplyr)
#   library(ggbio)
#   library(GenomicRanges)
#   library(Homo.sapiens) # is based on hg19
#   data("ideoCyto", package = "biovizBase")
#   data(genesymbol, package = "biovizBase")
#   library(BSgenome.Hsapiens.UCSC.hg19)
#   source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
#   
#   thisstart = max(0, thisstart)
#   thisend = min(thisend, seqlengths(BSgenome.Hsapiens.UCSC.hg19)[[thischr]])
#   
#   
#   ### BUILD TRANSLOCATION TRACK
#   
#   tx_original_2 = tx_original
#   colnames(tx_original_2)[4:7] = c("ChrB", "PosB", "ChrA", "PosA") 
#   tx_original_double = rbind(tx_original, tx_original_2)
#   tx_granges = tx_original_double %>% 
#     filter(Sample==thissample, SVCaller ==thissvcaller, ChrA == thischr, hasOverlap(PosA, PosA, thisstart, thisend)) %>% 
#     dplyr::select(ChrA, PosA, PosA, CircleGenomeClass)
#   tx_granges$Col = "bystander_tx"
#   tx_granges[tx_granges$Pos == txptinfo[[i, "TargetPos"]], "Col"] = "tx_of_interest"
#   if (nrow(tx_granges > 0)){
#     tx_granges = makeGRangesFromDataFrame(tx_granges, seqnames.field = "ChrA", start.field="PosA", end.field="PosA", ignore.strand = T, keep.extra.columns = T)
#     tx.track = autoplot(tx_granges, geom="rect", size=0.1, aes(color=Col), size=100) + scale_color_manual(values = c("bystander_tx" = "black", "tx_of_interest" = "red")) + guides(color=F) + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
#   } else {
#     tx_granges = GRanges()
#     tx.track = ggplot() + theme_null()
#   }
#   
#   
#   ### BUILD PALM TREE TRACK
#   
#   pt_granges = palmtrees %>% 
#     filter(Sample==thissample, SVCaller==thissvcaller) %>% 
#     filter(Sample==thissample, Chr == thischr, hasOverlap(FirstElement, LastElement, thisstart, thisend)) %>% 
#     dplyr::select(Chr, FirstElement, LastElement)
#   if (nrow(pt_granges) > 0){
#     pt_granges = makeGRangesFromDataFrame(pt_granges, seqnames.field = "Chr", start.field="FirstElement", end.field="LastElement", ignore.strand = T)
#     pt.track = autoplot(pt_granges, geom="rect", color=NA, fill="#E41A1C") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
#   } else {
#     pt_granges = GRanges()
#     pt.track = ggplot() + theme_null()    
#   }
#   
#   
#   ### BUILD CIRCLE-SEQ eccDNA
#   
#   circ_granges = circles %>%
#     filter(Method == "CircleSeq", CircleClass == "eccDNA")
#   filter(Sample==thissample, Chr == thischr, hasOverlap(Start, End, thisstart, thisend)) %>% 
#     dplyr::select(CircleChr, CircleStart, CircleEnd)
#   if (nrow(circ_granges)>0){
#     circ_granges = makeGRangesFromDataFrame(circ_granges, seqnames.field = "Chr", start.field="Start", end.field="End", ignore.strand = T)
#     circ.ecc.track = autoplot(circ_granges, geom="rect", color=NA, fill="#377EB8") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
#   } else {
#     circ.granges = GRanges()
#     circ.ecc.track = ggplot() + theme_null()
#   }
#   
#   ### BUILD CIRCLE-SEQ eccDNA
#   circ_granges = circles %>%
#     filter(Method == "CircleSeq", CircleClass == "ecDNA")
#   filter(Sample==thissample, Chr == thischr, hasOverlap(Start, End, thisstart, thisend)) %>% 
#     dplyr::select(CircleChr, CircleStart, CircleEnd)
#   if (nrow(circ_granges)>0){
#     circ_granges = makeGRangesFromDataFrame(circ_granges, seqnames.field = "Chr", start.field="Start", end.field="End", ignore.strand = T)
#     circ.ec.track = autoplot(circ_granges, geom="rect", color=NA, fill="#377EB8") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
#   } else {
#     circ.granges = GRanges()
#     circ.ec.track = ggplot() + theme_null()
#   }  
#   
#   # seqlevels(circ_granges, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
#   # seqinfo(circ_granges) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
#   # circ_granges <- trim(circ_granges)
#   
#   
#   ### BUILD WGS eccDNA TRACK
#   wgscirc_granges = merged_wgscircles %>%
#     filter(Sample==thissample, MergedCirclesChr == thischr, hasOverlap(MergedCirclesStart, MergedCirclesEnd, thisstart, thisend)) %>% 
#     dplyr::select(MergedCirclesChr, MergedCirclesStart, MergedCirclesEnd)
#   if (nrow(wgscirc_granges) > 0){
#     wgscirc_granges = makeGRangesFromDataFrame(wgscirc_granges, seqnames.field = "MergedCirclesChr", start.field="MergedCirclesStart", end.field="MergedCirclesEnd", ignore.strand = T)
#     wgscirc.track = autoplot(wgscirc_granges, geom="rect", color=NA, fill= "#4DAF4A") + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
#   }else{
#     wgscirc_granges = GRanges()
#     wgscirc.track = ggplot() + theme_null()
#   }
#   
#   
#   ### BUILD COPY NUMBER TRACK
#   
#   ascat_granges = ascat_cnv %>%
#     filter(Sample==thissample, Chr == thischr, hasOverlap(Start, End, thisstart, thisend)) %>% 
#     mutate(y = log10((TumorTotalCopyNumber + 0.01) / (NormalTotalCopyNumber + 0.01))) %>%
#     dplyr::select(Chr, Start, End, y)
#   if (nrow(ascat_granges)>0){
#     ascat_granges = makeGRangesFromDataFrame(ascat_granges, seqnames.field = "Chr", start.field="Start", end.field="End", ignore.strand = T, keep.extra.columns=T)
#     ascat.track = autoplot(ascat_granges, geom="bar", aes(y=y), fill="#984EA3", color="#984EA3", ylab="", size=0.1) + xlim(thisstart,thisend) + theme_clear() + theme(axis.line.x = element_blank())
#   } else {
#     ascat_granges = GRanges()
#     ascat.track = ggplot() + theme_null()
#   }
#   
#   
#   ### BUILD TRANSCRIPTS TRACK
#   
#   if (((thisend-thisstart)<=10000000) & ((thisend-thisstart)>20)){
#     gene.track = ggplot(data=data.frame(x = c(thisstart, thisend), y=c(0,0)), aes(x=x, y=y)) + geom_point(alpha=0) + xlim(thisstart, thisend) + annotation_custom(grid::textGrob("There is probably no \n transcripts in the region of interest"), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_clear() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) 
#     gene.track = try(autoplot(Homo.sapiens, 
#                               which = GRanges(seqnames=thischr,
#                                               ranges=IRanges(start = thisstart, end = thisend),
#                                               strand="*"),
#                               names.expr="SYMBOL",
#                               label.color="black") +
#                        xlim(thisstart,thisend) + 
#                        theme_clear(),    
#                      silent=TRUE)
#     if (is.character(gene.track) | isS4(gene.track)){ # if the try caught an error 
#       gene.track = ggplot(data=data.frame(x = c(thisstart, thisend), y=c(0,0)), aes(x=x, y=y)) + geom_point(alpha=0) + xlim(thisstart, thisend) + annotation_custom(grid::textGrob("There is probably no \n transcripts in the region of interest"), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_clear() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank()) 
#     }
#   } else {
#     gene.track=ggplot(data=data.frame(x = c(thisstart, thisend), y=c(0,0)), aes(x=x, y=y)) + geom_point(alpha=0) + xlim(thisstart, thisend) + annotation_custom(grid::textGrob("Too small or \n too large region \n to retrieve transcripts."), xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) + theme_clear() + theme(axis.line.x = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
#   }
#   
#   
#   ### PUT ALL TRACKS TOGETHER
#   
#   if (thissample %in% cscircles$Sample){
#     tracks("Translocations" = tx.track, "Palm Tree \n Region" = pt.track, "eccDNA \n (Circle-Seq)" = circ.track, "eccDNA \n (WGS)" = wgscirc.track, "Copy Number \n logFC" = ascat.track, "Transcripts" = gene.track, 
#            heights = c(1,1,1,1,1,2),
#            track.bg.color=NULL,
#            xlab=thischr,
#            title = paste0("Target ", txptinfo[i, "TargetID"]),
#            label.text.cex = 1,
#            label.bg.fill = "white",
#            label.bg.color = "white") +
#       theme(axis.line=element_line(color="black")) +
#       ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/StackFiguresTargetSites/", txptinfo[i,"TargetID"], "_StackFigure.pdf"), height=10, width=10)
#   } else {
#     tracks("Translocations" = tx.track, "Palm Tree \n Region" = pt.track, "eccDNA \n (WGS)" = wgscirc.track, "Copy Number \n logFC" = ascat.track, "Transcripts" = gene.track, 
#            heights = c(1,1,1,1,2),
#            track.bg.color=NULL,
#            xlab=thischr,
#            title = paste0("Target ", txptinfo[i, "TargetID"]),
#            label.text.cex = 1,
#            label.bg.fill = "white",
#            label.bg.color = "white") +
#       theme(axis.line=element_line(color="black")) +
#       ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/StackFiguresTargetSites/", txptinfo[i,"TargetID"], "_StackFigure.pdf"), height=10, width=10)
#   }
#   
# }
# 
# ## PLOT PALM TREE REGIONS
# for (i in i:nrow(palmtrees)){
#   palmtreelength = palmtrees[[i,"LastElement"]] - palmtrees[[i,"FirstElement"]]
#   plot_palm_tree(palmtrees[[i, "PalmTreeID"]],
#                  palmtrees[[i,"Sample"]], 
#                  palmtrees[[i,"SVCaller"]], 
#                  palmtrees[[i,"Chr"]], 
#                  palmtrees[[i,"FirstElement"]] - ceiling(0.25 * palmtreelength), 
#                  palmtrees[[i,"LastElement"]] +  ceiling(0.25 * palmtreelength), 
#                  palmtrees, tx_original, cscircles, merged_wgscircles, ascat_cnv, circlegenomeclass)
# }
# 
# ## PLOT TARGET REGIONS
# txptinfo$TargetID = paste0(txptinfo$PalmTreeID, "_to_", txptinfo$TargetChrom, ":", txptinfo$TargetPos)
# for (i in 8:nrow(txptinfo)){
#   try(plot_target_site(txptinfo, i, 
#                        txptinfo[i, "PalmTreeID"], 
#                        txptinfo[i, "Sample"], 
#                        txptinfo[i, "SVCaller"], 
#                        txptinfo[i, "TargetChrom"], 
#                        txptinfo[i, "TargetPos"]-500000, 
#                        txptinfo[i, "TargetPos"]+500000, 
#                        palmtrees, tx_original, cscircles, merged_wgscircles, ascat_cnv))
# }
# 
# 
# # library(plyranges)
# # this_sample = "NB2025"
# # this_svcaller = "Union"
# # 
# # this_txptinfo = txptinfo %>% filter(Sample == this_sample, SVCaller == this_svcaller)
# # 
# # this_tx_original = tx_original %>% filter(Sample == this_sample, SVCaller == this_svcaller)
# # 
# # this_tx_original=
# #   this_tx_original %>%
# #   mutate(
# #     isInPT =  
# #       (paste0(ChrA, ":", PosA) %in% paste0(this_txptinfo$PalmTreeChrom, ":", this_txptinfo$PalmTreePos)) |
# #       (paste0(ChrB, ":", PosB) %in% paste0(this_txptinfo$PalmTreeChrom, ":", this_txptinfo$PalmTreePos)))
# # 
# # 
# # this_tx_original_gr = makeGRangesFromDataFrame(
# #   this_tx_original,
# #   seqnames.field = "ChrA", start.field = "PosA", end.field = "PosA", 
# #   seqinfo = seqinfo(hg19sub), 
# #   keep.extra.columns = T)
# # LinkedTo = makeGRangesFromDataFrame(
# #   this_tx_original,
# #   seqnames.field = "ChrB", start.field = "PosB", end.field = "PosB", 
# #   seqinfo = seqinfo(hg19sub))
# # this_tx_original_gr$LinkedTo = LinkedTo
# # 
# # this_ascat = ascat_cnv %>%
# #   filter(Sample==this_sample) %>% 
# #   mutate(y = log2((TumorTotalCopyNumber + 0.01) / (NormalTotalCopyNumber + 0.01))) %>%
# #   dplyr::select(Chr, Start, End, y)
# # this_ascat_gr = makeGRangesFromDataFrame(this_ascat, seqnames.field = "Chr", start.field="Start", end.field="End", ignore.strand = T, keep.extra.columns=T, seqinfo = seqinfo(hg19sub))
# # 
# # bins <- tileGenome(seqinfo(hg19sub), tilewidth=500000, cut.last.tile.in.chrom=TRUE)
# # this_ascat_gr_binned = binnedAverage(bins=bins, coverage(this_ascat_gr, weight=this_ascat_gr$y), "y")
# # 
# # 
# # this_palmtrees_gr = palmtrees_gr %>% filter(Sample == this_sample, SVCaller == this_svcaller)
# # seqinfo(this_palmtrees_gr) = seqinfo(keepStandardChromosomes(GRangesForBSGenome("hg19"), pruning.mode = "coarse"))
# # 
# # 
# # ggbio() + 
# #   circle(this_tx_original_gr, geom ="link", linked.to="LinkedTo", aes(color=isInPT), size=0.3) +
# #   circle(this_palmtrees_gr, geom = "rect", fill = "red", color="red") +
# #   circle(this_ascat_gr_binned, geom="bar", aes(y=y, fill=y), color=NA, trackWidth=6, grid=T, 
# #          grid.background=NA, grid.line = "gray90", grid.n=3) +
# #   scale_fill_gradient2(low="blue", high="red", mid="gray90") + 
# #   #circle(hg19sub, geom = "ideo", fill = "white") +
# #   circle(hg19sub, geom = "scale", size = 2) +
# #   circle(hg19sub, geom = "text", aes(label = seqnames), size = 5, vjust=-1) +
# #   scale_color_manual(values = c("grey80", "red"))
# # 
# # hg19sub = keepStandardChromosomes(GRangesForBSGenome("hg19"), pruning.mode = "coarse")
# # 
# # 
# # bins <- tileGenome(seqinfo(hg19sub), tilewidth=10000000, cut.last.tile.in.chrom=TRUE)
# # this_ascat_gr_binned = binnedAverage(bins=bins, coverage(this_ascat_gr, weight=this_ascat_gr$y), "y")
