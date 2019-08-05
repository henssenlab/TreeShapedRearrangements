rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
source("~/Desktop/PalmTrees/Analysis/Code/ParseSVCallerData.R")
source("~/Desktop/PalmTrees/Analysis/Code/AlmostDuplicateBreakpoints.R")
source("~/Desktop/PalmTrees/Analysis/Code/CallHighDensityRegion.R")

tx = parse_pedpancan("~/Desktop/PalmTrees/Data/PedPanCanSVs.csv")
meta = read.table("~/Desktop/PalmTrees/Data/PedPanCanMeta.csv", sep=";", header=T)
wgs_samples = meta %>% filter(seq_type == "wgs") %>% .$sample
tx = tx %>% filter(Sample %in% wgs_samples)

# tx %>% dplyr::select(Sample, ChrA, PosA, ChrB, PosB) %>% write.table("~/Desktop/PalmTrees/Analysis/Shiny/PalmTrees/pedpancan_wgs_txcalls.tsv",
#                                                                      sep = "\t",
#                                                                      quote=F,
#                                                                      row.names=F,
#                                                                      col.names=T)

#######################  MANUAL FILTERING  ######################                  

# Where are the Translocations across all Samples? Hotspots?
data.frame(Chr = c(tx$ChrA, tx$ChrB), Pos = c(tx$PosA, tx$PosB)) %>%
  ggplot(aes(x=Pos)) +
  facet_wrap(Chr~.) + 
  geom_histogram(bins=100000)

# Translocations per Sample
library(ggrepel)
tx %>%
  group_by(Sample) %>%
  summarise(ntx=n()) %>% 
  ungroup() %>%
  mutate(Rank = row_number(ntx)) %>%
  ggplot(aes(x = ntx, y=Rank)) +
  geom_point() 

zscore = function(v) ( (v-mean(v, na.rm=T)) / (sd(v, na.rm=T)+0.00001) )

# OutlierSamples = tx %>%
#   group_by(Sample) %>%
#   summarise(ntx=n()) %>%
#   ungroup() %>%
#   mutate(z = zscore(ntx)) %>%
#   filter(z >= 2) %>%
#   .$Sample
# 
# tx = tx %>% filter(!(Sample %in% OutlierSamples))

tx = tx %>% filter(Sample != "sjos103")

tx$Sample  = droplevels(tx$Sample)

#################################################################

#### BEGIN PALM TREE CALLING

# keep the original, unfiltered translocation data as tx_original
tx_original = tx

# filter chromosome pairs with many tx
'%notin%' <- function(x,y)!('%in%'(x,y))
tx$ChrTuple = as.character(
  mapply(function (a,b) do.call(paste0, as.list(sort(c(a, b)))),
         tx$ChrA, tx$ChrB)
)
tx$SampleSVCallerChrTuple = paste0(tx$Sample, tx$SVCaller, tx$ChrTuple)
filteroutsamplechrtuples = tx %>% group_by(SampleSVCallerChrTuple) %>% summarise(n = n()) %>% filter(n>=5) %>% .$SampleSVCallerChrTuple
tx = tx %>% filter(SampleSVCallerChrTuple %notin% filteroutsamplechrtuples)

# txdouble binds together tx plus tx with Breakpoints A and B swapped  
tx2 = tx
colnames(tx2)[2:5] = c("ChrB", "PosB", "ChrA", "PosA") 
txdouble_unfiltered = rbind(tx, tx2)
rm(tx2)

# take interchromosomal tx only
tx = tx[tx$ChrA != tx$ChrB,]
tx2 = tx
colnames(tx2)[2:5] = c("ChrB", "PosB", "ChrA", "PosA") 
txdouble = rbind(tx, tx2) # That contains a row for each breakpoint instead of a row for each breakpoint pair
txdouble$SampleSVCallerChrA = paste0(txdouble$Sample, ":", txdouble$SVCaller, ":", txdouble$ChrA)
txdouble$SampleSVCaller = paste0(txdouble$Sample, ":", txdouble$SVCaller)
length_txdouble_first = nrow(txdouble)
rm(tx2)

# filter out similar breakpoint pairs
filteredtxdouble = txdouble[0,]
samplesXsvcallers = unique(txdouble$SampleSVCaller)
for (i in 1:length(samplesXsvcallers)){
  this_samplexsvcaller_data = txdouble %>% filter(SampleSVCaller == samplesXsvcallers[i])
  chrs = unique(c(this_samplexsvcaller_data$ChrA, this_samplexsvcaller_data$ChrB))
  for (chrA_i in 1:length(chrs)){
    for (chrB_i in 1:length(chrs)){
      this_chr_data = this_samplexsvcaller_data %>% filter(ChrA == chrs[chrA_i], ChrB == chrs[chrB_i])
      aldups = almost_duplicate_breakpoints(as.vector(this_chr_data$PosA), as.vector(this_chr_data$PosB), 10000000)
      if (length(aldups) > 0) this_chr_data = this_chr_data[-almost_duplicate_breakpoints(as.vector(this_chr_data$PosA), as.vector(this_chr_data$PosB), 10000000),]
      filteredtxdouble = rbind(filteredtxdouble, this_chr_data)
    }
  }
}
txdouble = filteredtxdouble
rm(this_chr_data, this_samplexsvcaller_data, chrs, aldups, filteredtxdouble)

# filter out shattered chromosomes (>50 translocations on one chromosome)
shattering_threshold = 25
samplesXsvcallers = unique(txdouble$SampleSVCaller)
for (i in 1:length(samplesXsvcallers)){
  this_samplexsvcaller_data = txdouble %>% filter(SampleSVCaller == samplesXsvcallers[i])
  chrs = unique(this_samplexsvcaller_data$ChrA)
  for (chr_i in 1:length(chrs)){
    if (nrow(this_samplexsvcaller_data %>% filter(ChrA == chrs[chr_i] | ChrB == chrs[chr_i])) > 2*shattering_threshold){
      txdouble = txdouble %>% filter(!(SampleSVCaller == samplesXsvcallers[i]) | !(ChrA == chrs[chr_i] | ChrB == chrs[chr_i]))
    }
  }
}
rm(this_samplexsvcaller_data)
length_txdouble_afterchromosomeshattering = nrow(txdouble)

windowlength = 4000000
resol = 25000
threshold = 3

palmtrees = data.frame(IntervalID=character(), Cohort=character(), SVCaller=character, Sample=character(), Chr=character(), Start=double(), End=double(), PalmTreeID=character())
samples = unique(txdouble$Sample)
svcallers = unique(txdouble$SVCaller)
for (sample_i in 1:length(samples)){
  print(samples[sample_i])
  for (svcaller_i in 1:length(svcallers)){
    d1 = txdouble %>% filter(Sample == samples[sample_i], SVCaller == svcallers[svcaller_i])
    chrs = unique(d1$ChrA)
    for (chr_i in 1:length(chrs)){
      d = d1 %>% filter(ChrA == chrs[chr_i])
      d = as.data.frame(d$PosA)
      colnames(d) = c("PosA")
      #d = d %>% distinct() # KH 18 12 05
      thispts = call_high_density_region(d$PosA, windowlength, resol, threshold)
      if (length(thispts) > 0){
        thispts$Sample = samples[sample_i]
        thispts$SVCaller = svcallers[svcaller_i]
        thispts$Chr = chrs[chr_i]
        thispts$PalmTreeID = paste0(thispts$Sample, "_", thispts$SVCaller, "_", thispts$Chr, ":", as.character(thispts$FirstElement), "-", as.character(thispts$LastElement))
        palmtrees = rbind(palmtrees, thispts)
      }
    }
  }
}

# exclude Palmtrees e.g. on mitochondria or undefined assembly sequences
palmtrees = palmtrees %>% filter(isdefchrom(Chr))
  
txdouble = txdouble_unfiltered
txdouble$PalmTreeID = NA
txdouble_ptinfo = txdouble[0,]
txptinfo = tx_original[0,]
txdouble$isPalmTree = FALSE
for (i in 1:nrow(palmtrees)){
  txptinfo = rbind(txptinfo, 
                   tx_original %>% 
                     filter(Sample == palmtrees[[i, "Sample"]], SVCaller == palmtrees[[i, "SVCaller"]],
                            ChrA == palmtrees[[i, "Chr"]], PosA >= palmtrees[[i, "FirstElement"]], PosA <= palmtrees[[i, "LastElement"]]) %>%
                     mutate(PalmTreeChrom=ChrA, PalmTreePos=PosA, TargetChrom=ChrB, TargetPos=PosB, PalmTreeID = palmtrees[[i, "PalmTreeID"]]),
                   tx_original %>% 
                     filter(Sample == palmtrees[[i, "Sample"]], SVCaller == palmtrees[[i, "SVCaller"]], 
                            ChrB == palmtrees[[i, "Chr"]], PosB >= palmtrees[[i, "FirstElement"]], PosB <= palmtrees[[i, "LastElement"]], 
                            ((ChrA != palmtrees[[i, "Chr"]]) | (PosA < palmtrees[[i, "FirstElement"]]) | (PosA > palmtrees[[i, "LastElement"]]))) %>% 
                     mutate(PalmTreeChrom=ChrB, PalmTreePos=PosB, TargetChrom=ChrA, TargetPos=PosA, PalmTreeID = palmtrees[[i, "PalmTreeID"]]))
  
  txdouble[(txdouble$Sample == palmtrees[[i, "Sample"]]) & (txdouble$SVCaller == palmtrees[[i, "SVCaller"]]) & (txdouble$ChrA == palmtrees[[i, "Chr"]]) & (txdouble$PosA >= palmtrees[[i, "FirstElement"]]) & (txdouble$PosA <= palmtrees[[i, "LastElement"]]), "isPalmTree"] = TRUE
  txdouble[(txdouble$Sample == palmtrees[[i, "Sample"]]) & (txdouble$SVCaller == palmtrees[[i, "SVCaller"]]) & (txdouble$ChrB == palmtrees[[i, "Chr"]]) & (txdouble$PosB >= palmtrees[[i, "FirstElement"]]) & (txdouble$PosB <= palmtrees[[i, "LastElement"]]), "isPalmTree"] = TRUE
}

txptinfo = txptinfo %>% dplyr::select(-ChrA, -PosA, -ChrB, -PosB)

palmtrees = palmtrees %>% dplyr::select(-Start, -End, -StartContribs, -EndContribs)

# make GRanges Object for further use
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
palmtrees_gr = makeGRangesFromDataFrame(as_tibble(palmtrees),
                                        seqnames.field = "Chr",
                                        start.field = "FirstElement",
                                        end.field = "LastElement",
                                        ignore.strand = T,
                                        keep.extra.columns = T)
seqlevels(palmtrees_gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)[1:23]
save.image(file = "~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan.Rdata")

# load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan.Rdata")
# 
# library(circlize)
# 
# plot_onesample_onesvcaller = function(palmtrees, txdouble_ptinfo, txdouble, sample, svcaller, fname){
#   
#   mycols = rand_color(n=100, luminosity = "dark")
#   
#   palmtrees_of_interest = palmtrees %>% filter(Sample == sample, SVCaller == svcaller)
#   
#   pdf(file=fname)
#   circos.clear()
#   circos.par("start.degree" = 90)
#   #circos.initializeWithIdeogram(species = "hg19")
#   circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
#   text(0, 0, "", cex = 1) 
#   
#   palmtrees_bed = palmtrees_of_interest[,c("Chr", "FirstElement", "LastElement")]
#   
#   # this is just for plotting
#   diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
#   if (sum(diff<5000000)>0){
#     palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff)/2
#     palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff)/2
#   }
#   
#   if (nrow(palmtrees_bed)>0){
#     palmtrees_bed$value1 = 1:nrow(palmtrees_of_interest)
#     palmtrees_bed$value2 = palmtrees_of_interest$PalmTreeID
#   }else{
#     palmtrees_bed$value1 = data.frame()
#     palmtrees_bed$value2 = data.frame()
#   }
#   colnames(palmtrees_bed) = c("chr", "start", "end", "value1", "value2")
#   
#   #bed_palmtree_A = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
#   #bed_palmtree_B = data.frame("ChrA"=c(), "PosA"=c(), "PosA"=c(), "PalmTreeIndex"=c())
#   
#   bed_palmtree_A = data.frame()
#   bed_palmtree_B = data.frame()
#   if (nrow(palmtrees_bed) > 0){
#     for (i in 1:nrow(palmtrees_bed)){
#       palmtree_tx = txdouble_ptinfo %>% filter(PalmTreeID == as.character(palmtrees_bed[i, "value2"]))
#       palmtree_tx$PalmTreeIndex = i
#       palmtree_tx = palmtree_tx %>% dplyr::select(PalmTreeChrom, PalmTreePos, TargetChrom, TargetPos, PalmTreeIndex) %>% distinct()
#       temp = palmtree_tx[,c("PalmTreeChrom", "PalmTreePos", "PalmTreePos", "PalmTreeIndex")]
#       colnames(temp) = c("chr", "start", "end", "value1")
#       bed_palmtree_A = rbind(bed_palmtree_A, temp)
#       colnames(bed_palmtree_A) = c("chr", "start", "end", "value1")
#       temp = palmtree_tx[,c("TargetChrom", "TargetPos", "TargetPos", "PalmTreeIndex")]
#       colnames(temp) = c("chr", "start", "end", "value1")
#       bed_palmtree_B = rbind(bed_palmtree_B, temp)
#       colnames(bed_palmtree_B) = c("chr", "start", "end", "value1")
#     }
#   }
#   
#   circos.genomicTrack(as.data.frame(palmtrees_bed),
#                       panel.fun = function(region, value, ...){
#                         circos.genomicRect(region, value, col=mycols[value[[1]]], border=NA, ...)
#                       },
#                       ylim=c(0,0.2), track.height=0.33*circos.par("track.height"))
#   
#   nopalmtree_tx = tx_original %>% filter(Sample == sample, SVCaller == svcaller, isdefchrom(ChrA), isdefchrom(ChrB))
#   bed_nopalmtree_A = nopalmtree_tx[,c("ChrA", "PosA", "PosA")]
#   if (nrow(bed_nopalmtree_A) > 0){
#     bed_nopalmtree_A$value1 = 0
#     colnames(bed_nopalmtree_A) = c("chr", "start", "end", "value1")
#   }
#   
#   bed_nopalmtree_B = nopalmtree_tx[,c("ChrB", "PosB", "PosB")]
#   if (nrow(bed_nopalmtree_B) > 0){
#     bed_nopalmtree_B$value1 = 0
#     colnames(bed_nopalmtree_B) = c("chr", "start", "end", "value1")
#   }
#   circos.genomicLink(bed_nopalmtree_A, bed_nopalmtree_B, col="lightgray")
#   
#   if (nrow(palmtrees_bed) > 0) circos.genomicLink(bed_palmtree_A, bed_palmtree_B, col=mycols[bed_palmtree_A[,"value1"]])
#   
#   title(paste(sample, svcaller))
#   dev.off()
#   circos.clear()
# }
# 
# samples = unique(as.character(tx$Sample))
# svcallers = unique(as.character(tx$SVCaller))
# for (sai in 1:length(samples)){
#   print(samples[sai])
#   for (svi in 1:length(svcallers)){
#     if (nrow(tx %>% filter(Sample == samples[sai], SVCaller == svcallers[svi]))>0){
#       plot_onesample_onesvcaller(palmtrees, txptinfo, txdouble, samples[sai], svcallers[svi], paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlotsPedPanCan/", as.character(svcallers[svi]), "/", as.character(samples[sai]), as.character(svcallers[svi]), "PalmTreeCircos.pdf"))
#     }
#   }
# }
# 
