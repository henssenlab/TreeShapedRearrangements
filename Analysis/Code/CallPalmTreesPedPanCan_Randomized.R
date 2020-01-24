rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
source("~/Desktop/PalmTrees/Analysis/Code/ParseSVCallerData.R")
source("~/Desktop/PalmTrees/Analysis/Code/AlmostDuplicateBreakpoints.R")
source("~/Desktop/PalmTrees/Analysis/Code/CallHighDensityRegion.R")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")
library(BSgenome.Hsapiens.UCSC.hg19)

tx = parse_pedpancan("~/Desktop/PalmTrees/Data/PedPanCanSVs.csv")
meta = read.table("~/Desktop/PalmTrees/Data/PedPanCanMeta.csv", sep=";", header=T)
wgs_samples = meta %>% filter(seq_type == "wgs") %>% .$sample
tx = tx %>% filter(Sample %in% wgs_samples)
tx = tx %>% filter(Sample != "sjos103")

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
#rm(tx2)

# take interchromosomal tx only
tx = tx[tx$ChrA != tx$ChrB,]
tx2 = tx
colnames(tx2)[2:5] = c("ChrB", "PosB", "ChrA", "PosA") 
txdouble = rbind(tx, tx2) # That contains a row for each breakpoint instead of a row for each breakpoint pair
txdouble$SampleSVCallerChrA = paste0(txdouble$Sample, ":", txdouble$SVCaller, ":", txdouble$ChrA)
txdouble$SampleSVCaller = paste0(txdouble$Sample, ":", txdouble$SVCaller)
length_txdouble_first = nrow(txdouble)
#rm(tx2)

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
#rm(this_chr_data, this_samplexsvcaller_data, chrs, aldups, filteredtxdouble)

# filter out shattered chromosomes (>50 translocations on one chromosome) --------------------------- this does nothing in PedPanCancer
shattering_threshold = 25
samplesXsvcallers = unique(txdouble$SampleSVCaller)
filtered_samplechrs = list()
k = 1
for (i in 1:length(samplesXsvcallers)){
  this_samplexsvcaller_data = txdouble %>% filter(SampleSVCaller == samplesXsvcallers[i])
  chrs = unique(this_samplexsvcaller_data$ChrA)
  for (chr_i in 1:length(chrs)){
    if (nrow(this_samplexsvcaller_data %>% filter(ChrA == chrs[chr_i] | ChrB == chrs[chr_i])) > 2*shattering_threshold){
      print(chrs[chr_i])
      txdouble = txdouble %>% filter(!(SampleSVCaller == samplesXsvcallers[i]) | !(ChrA == chrs[chr_i] | ChrB == chrs[chr_i]))
      filtered_samplechrs[[k]] = paste0(samplesXsvcallers[i], ":", chrs[chr_i])
      k=k+1
    }
  }
}
#rm(this_samplexsvcaller_data)
length_txdouble_afterchromosomeshattering = nrow(txdouble)
filtered_samplechrs = do.call(c, filtered_samplechrs)

old_txdouble = txdouble

## Randomize
new_tx_all_samples = list()
for (i in 1:length(samplesXsvcallers)){
  this_samplexsvcaller_data = txdouble %>% filter(SampleSVCaller == samplesXsvcallers[i])
  if ((nrow(this_samplexsvcaller_data) %% 2) != 0) stop("this_samplexsvcaller_data does not have even number of rows")
  this_samplexsvcaller_ntx = nrow(this_samplexsvcaller_data) / 2

  new_tx = data.frame()
  while (nrow(new_tx) < this_samplexsvcaller_ntx){
    A = sample_random_regions_broad(50*this_samplexsvcaller_ntx, 0) %>% dplyr::select(Chr, Start) %>% dplyr::rename(ChrA = Chr, PosA = Start)
    B = sample_random_regions_broad(50*this_samplexsvcaller_ntx, 0) %>% dplyr::select(Chr, Start) %>% dplyr::rename(ChrB = Chr, PosB = Start)
    
    new_tx = bind_cols(A,B) %>%
      filter(ChrA != ChrB)
    
    new_tx$ChrA = as.character(new_tx$ChrA)
    new_tx$PosA = new_tx$PosA
    new_tx$ChrB = as.character(new_tx$ChrB)
    new_tx$PosB = new_tx$PosB

    new_tx$SampleChrA = paste0(this_samplexsvcaller_data[[1,"SampleSVCaller"]],":",new_tx$ChrA)
    new_tx$SampleChrB = paste0(this_samplexsvcaller_data[[1,"SampleSVCaller"]],":",new_tx$ChrB)
    
    new_tx$ChrTuple = as.character(
      mapply(function (a,b) do.call(paste0, as.list(sort(c(a, b)))),
             new_tx$ChrA, new_tx$ChrB)
    )
    new_tx$Sample = this_samplexsvcaller_data[[1,"Sample"]]
    new_tx$SVCaller = this_samplexsvcaller_data[[1,"SVCaller"]]
    new_tx$SampleChrTuple = paste0(new_tx$Sample, new_tx$SVCaller,new_tx$ChrTuple)
    
    new_tx = new_tx %>% filter(SampleChrTuple %notin% filteroutsamplechrtuples,
                               SampleChrA %notin% filtered_samplechrs,
                               SampleChrB %notin% filtered_samplechrs)
  }
  new_tx = new_tx[1:this_samplexsvcaller_ntx,]
  new_tx$ChrA = as.character(new_tx$ChrA)
  new_tx$PosA = new_tx$PosA
  new_tx$ChrB = as.character(new_tx$ChrB)
  new_tx$PosB = new_tx$PosB

  new_tx_all_samples[[i]] = new_tx
}
new_tx_all_samples = do.call(rbind, new_tx_all_samples) 
new_tx_all_samples = new_tx_all_samples %>% dplyr::select(-SampleChrA, -SampleChrB)

tx = new_tx_all_samples
tx_original = new_tx_all_samples

new_tx_all_samples_2 = new_tx_all_samples
colnames(new_tx_all_samples_2)[1:4] = c("ChrB", "PosB", "ChrA", "PosA") 
txdouble = rbind(new_tx_all_samples, new_tx_all_samples_2)

# nrow(txdouble)
# nrow(old_txdouble)
# 
# old_txdouble$Sample %>% table()
# txdouble$Sample %>% table()
# 
# old_txdouble %>%
#   ggplot(aes(x=PosA)) +
#   geom_density()+
#   facet_wrap(ChrA~.)
# 
# txdouble %>%
#   ggplot(aes(x=PosA)) +
#   geom_density()+
#   facet_wrap(ChrA~.)
# # pretty uniform distribution, is different
# 
# test_sample = "sjos013"
# 
# txdouble %>% filter(Sample == test_sample) %>% View()
# 
# filteroutsamplechrtuples[grepl(test_sample, filteroutsamplechrtuples)]
# txdouble %>% filter(Sample == test_sample) %>% .$ChrTuple %>% table()  # the forbidden tuples do not appear
# old_txdouble %>% filter(Sample == test_sample) %>% .$ChrTuple %>% table() # the distribution is different from before randomization
# txdouble %>% filter(Sample == test_sample) %>% nrow()
# old_txdouble %>% filter(Sample == test_sample) %>% nrow() # it is the same number of rearrangements
# 
# old_txdouble %>% filter(Sample == test_sample) %>%
#   ggplot(aes(x=PosA)) +
#   geom_density()+
#   facet_wrap(ChrA~.)
# # we still see spikes
# 
# txdouble %>% filter(Sample == test_sample) %>%
#   ggplot(aes(x=PosA)) +
#   geom_density()+
#   facet_wrap(ChrA~.)
# # pretty uniform distribution


########## SLIDING WINDOW
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
save.image(file = "~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan_RandomizedBreakpoints_Method3.Rdata")

## Explore randomization visually
rm(list=ls())   
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan.Rdata")
meta = meta %>% dplyr::rename(Sample = sample)
real_tx_original = tx_original
real_palmtrees = palmtrees
real_npalmtrees = real_palmtrees %>%
  left_join(meta) %>% 
  group_by(Sample, entity) %>% 
  summarise(npt_real = n_distinct(PalmTreeID)) %>% 
  ungroup() 

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan_RandomizedBreakpoints_Method3.Rdata")
meta = meta %>% dplyr::rename(Sample = sample)
npalmtrees = palmtrees %>%
  left_join(meta) %>% 
  group_by(Sample, entity) %>% 
  summarise(npt = n_distinct(PalmTreeID)) %>% 
  ungroup() 

compare_real_vs_permuted = full_join(real_npalmtrees, npalmtrees) %>%
  mutate(npt_real = ifelse(is.na(npt_real), 0, npt_real),
         npt = ifelse(is.na(npt), 0, npt)) %>%
  mutate(FC = (npt_real + 1) / (npt + 1)) %>% 
  mutate(hasRealPT = npt_real>0,
         hasSimulPT = npt>0) %>% 
  mutate(entity = as.character(entity)) %>% 
  mutate(entity = ifelse(grepl("mb_", entity), "mb", entity)) %>%
  mutate(entity = ifelse(grepl("hgg", entity), "hgg", entity)) %>%
  mutate(entity = ifelse(grepl("b-all", entity), "b-all", entity)) %>%
  mutate(entity = ifelse(grepl("epd_", entity), "epd", entity)) %>%
  mutate(entity = toupper(entity))

# Probability that a Sample has a palm tree in simulation if it has a palm tree
compare_real_vs_permuted %>% 
  summarise(sum_hasSimulPT = sum(hasSimulPT),
            sum_hasRealPT = sum(hasRealPT),
            sum_hasSimulPT_divby546 = sum_hasSimulPT/546,
            sum_hasRealPT_divby546 = sum_hasRealPT/546)

# by entity (Email Table 2)  <--------------------------------------------------------
compare_real_vs_permuted %>% 
  group_by(entity) %>% 
  summarise(sum_hasRealPT = sum(hasRealPT),
            sum_hasSimulPT = sum(hasSimulPT))

# a crude, FDR estimate
compare_real_vs_permuted %>% 
  filter(hasRealPT) %>% 
  summarise(mean(hasSimulPT))

# by entity (Email Table 1) <--------------------------------------------------------
compare_real_vs_permuted %>% 
  filter(hasRealPT) %>% 
  group_by(entity) %>% 
  summarise(prob = mean(hasSimulPT))


compare_real_vs_permuted %>% 
  group_by(entity) %>% 
  summarise(sum_npt = sum(npt),
            sum_npt_real = sum(npt_real))

#########

all_samples = 
  meta %>% filter(seq_type == "wgs") %>% dplyr::rename(Sample = sample) %>% 
  dplyr::select(Sample, entity) %>% distinct() %>% 
  left_join(npalmtrees) %>%
  mutate(npt = ifelse(is.na(npt), 0, npt)) %>%
  mutate(entity = as.character(entity)) %>% 
  mutate(entity = ifelse(grepl("mb_", entity), "mb", entity)) %>%
  mutate(entity = ifelse(grepl("hgg", entity), "hgg", entity)) %>%
  mutate(entity = ifelse(grepl("b-all", entity), "b-all", entity)) %>%
  mutate(entity = ifelse(grepl("epd_", entity), "epd", entity)) %>%
  mutate(entity = toupper(entity))

all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0)) %>% 
  gather("Class", "NSamples", 3:4) %>% 
  mutate(entity = forcats::fct_reorder(entity, -nsamples)) %>% 
  filter(nsamples > 0) %>% 
  mutate(entity = droplevels(entity)) %>% 
  mutate(Class = factor(Class, levels=c("nsamples_without_pt", "nsamples_with_pt"))) %>% 
  ggplot(aes(x=entity, y=NSamples, color=Class, fill=Class)) + 
  geom_col() + 
  theme_kons1() + 
  coord_flip()

all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0)) %>% 
  mutate(entity = forcats::fct_reorder(entity, percsamples_with_pt)) %>% 
  filter(nsamples > 0) %>% 
  mutate(entity = droplevels(entity)) %>% 
  ggplot(aes(x=entity, y=percsamples_with_pt)) + 
  geom_col() + 
  theme_kons1() + 
  coord_flip()


all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0))  

# load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan_RandomizedBreakpoints.Rdata")
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

