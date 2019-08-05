rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)

source("~/Desktop/PalmTrees/Analysis/Code/AlmostDuplicateBreakpoints.R")
source("~/Desktop/PalmTrees/Analysis/Code/CallHighDensityRegion.R")

#source("~/Desktop/PalmTrees/Analysis/Code/ParseMergedTx.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/MergedTx.Rdata")

tx %>% 
  filter(Cohort == "Peifer") %>%
  dplyr::select(SVCaller, Sample, ChrA, PosA, ChrB, PosB) %>%
  distinct() %>% 
  filter(Sample != "NBL50", Sample != "NBL49", Sample != "NBL47") %>% 
  write.table("~/Desktop/PalmTrees/Analysis/Shiny/PalmTrees/peifer_tx.tsv", sep="\t", quote=F, row.names=F, col.names=T)

tx %>% 
  filter(Cohort == "Berlin") %>%
  dplyr::select(SVCaller, Sample, ChrA, PosA, ChrB, PosB) %>%
  distinct() %>% 
  write.table("~/Desktop/PalmTrees/Analysis/Shiny/PalmTrees/berlin_tx.tsv", sep="\t", quote=F, row.names=F, col.names=T)

# #######################  MANUAL FILTERING ######################################                   
# # List by number of translocations
# tx %>% filter(SVCaller != "AtLeastTwo") %>% group_by(Sample, SVCaller) %>% summarise(n=n()) %>% arrange(desc(n))
# # Mean+SD of number of tx per sample 
# tx %>% filter(SVCaller != "AtLeastTwo") %>% group_by(Sample, SVCaller) %>% summarise(n=n()) %>% ungroup() %>% summarise(mean=mean(n), sd=sd(n))
# # Mean+SD of number of tx per sample for each caller
# tx %>% filter(SVCaller != "AtLeastTwo") %>% group_by(Sample, SVCaller) %>% summarise(n=n()) %>% ungroup() %>% group_by(SVCaller) %>% summarise(mean=mean(n), sd=sd(n))
# 
# # Exclude samples with exceptionally many translocations
# #tx = tx %>% filter(Sample != "NBL50", Sample != "NBL55", Sample != "NBL49", Sample != "NBL53", Sample != "NBL54", Sample != "NBL39", Sample != "NBL61")
# tx = tx %>% filter(Sample != "NBL50", Sample != "NBL49", Sample != "NBL47")
# #atleasttwo = atleasttwo %>% filter(Sample != "NBL50", Sample != "NBL49", Sample != "NBL47")
# 
# # comb = data.frame("SVCaller" = c(tx$SVCaller, tx$SVCaller),
# #                   "Chr" = c(tx$ChrA, tx$ChrB),
# #                   "Pos" = c(tx$PosA, tx$PosB))
# # comb$PosBin = paste0(comb$Chr, ":", ceiling(comb$Pos / 1000000))
# # comb %>% filter(SVCaller == "Union") %>% group_by(PosBin) %>% summarise(n=n()) %>% ungroup() %>% View()
# 
# # Exclude the weird locus around LINC00486
# tx = tx %>% filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
# #atleasttwo = atleasttwo %>% filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
# 
# # Exclude weird locus on chr12
# #tx = tx %>% filter(!(ChrA == "chr12" & PosA >= 66451367 & PosA <= 66451465), !(ChrB == "chr12" & PosB >= 66451367 & PosB <= 66451465))
# #atleasttwo = atleasttwo %>% filter(!(ChrA == "chr12" & PosA >= 66451367 & PosA <= 66451465), !(ChrB == "chr12" & PosB >= 66451367 & PosB <= 66451465))
# 
# # How many tx are called by different callers?
# tx %>% group_by(SVCaller) %>% summarise(n=n_distinct(BPID))

#### BEGIN PALM TREE CALLING

# keep the original, unfiltered translocation data as tx_original
tx_original = tx

# Exclude shattered chromosomes 3 and 5 in Sample NB2028 by hand
#tx = tx %>% filter((Sample != "NB2028") | !(ChrA %in% c("chr3", "chr5")) | !(ChrB %in% c("chr3", "chr5")))

# filter chromosome pairs with many tx
'%notin%' <- function(x,y)!('%in%'(x,y))
tx$ChrTuple = as.character(
  mapply(function (a,b) do.call(paste0, as.list(sort(c(a, b)))),
         tx$ChrA, tx$ChrB)
)
tx$SampleSVCallerChrTuple = paste0(tx$Sample, tx$SVCaller, tx$ChrTuple)
filteroutsamplesvcallerchrtuples = tx %>% group_by(SampleSVCallerChrTuple) %>% summarise(n = n_distinct(BPID)) %>% filter(n>=5) %>% .$SampleSVCallerChrTuple
tx = tx %>% filter(SampleSVCallerChrTuple %notin% filteroutsamplesvcallerchrtuples)

# txdouble binds together tx plus tx with Breakpoints A and B swapped  
tx2 = tx %>%
  mutate(
    ChrANew = ChrB,
    PosANew = PosB,
    ChrBNew = ChrA,
    PosBNew = PosA,
    ChrA = ChrANew,
    PosA = PosANew,
    ChrB = ChrBNew,
    PosB = PosBNew
  ) %>%
  dplyr::select(-ChrANew, -PosANew, -ChrBNew, -PosBNew)
txdouble_unfiltered = bind_rows(tx, tx2)
rm(tx2)

# take interchromosomal tx only
tx2 = tx %>%
  mutate(
    ChrANew = ChrB,
    PosANew = PosB,
    ChrBNew = ChrA,
    PosBNew = PosA,
    ChrA = ChrANew,
    PosA = PosANew,
    ChrB = ChrBNew,
    PosB = PosBNew
  ) %>%
  dplyr::select(-ChrANew, -PosANew, -ChrBNew, -PosBNew)
txdouble = bind_rows(tx, tx2)
# That contains a row for each breakpoint instead of a row for each breakpoint pair
txdouble$SampleSVCallerChrA = paste0(txdouble$Sample, ":", txdouble$SVCaller, ":", txdouble$ChrA)
txdouble$SampleSVCaller = paste0(txdouble$Sample, ":", txdouble$SVCaller)
length_txdouble_first = nrow(txdouble)
rm(tx2)

# filter out similar breakpoints
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

# exclude Chromosomes with more than 15 palm trees for one sv calling approach
# samples = unique(palmtrees$Sample)
# svcallers = unique(palmtrees$SVCaller)
# for (sample_i in 1:length(samples)){
#   for (svcaller_i in 1:length(svcallers)){
#     if (palmtrees %>% 
#         filter(Sample == samples[sample_i], SVCaller == svcallers[svcaller_i]) %>%
#         nrow() > 15){
#       palmtrees = palmtrees %>% filter(Sample != samples[sample_i] | SVCaller != svcallers[svcaller_i])
#     }
#   }
# }
  
txdouble = txdouble_unfiltered
txdouble$PalmTreeID = NA
txdouble_ptinfo = txdouble[0,]
txptinfo = tx_original[0,]
txdouble$isPalmTree = FALSE
for (i in 1:nrow(palmtrees)){
  txptinfo = rbind(txptinfo, 
                   tx_original %>% filter(Sample == palmtrees[[i, "Sample"]], SVCaller == palmtrees[[i, "SVCaller"]], ChrA == palmtrees[[i, "Chr"]], PosA >= palmtrees[[i, "FirstElement"]], PosA <= palmtrees[[i, "LastElement"]]) %>% mutate(PalmTreeChrom=ChrA, PalmTreePos=PosA, TargetChrom=ChrB, TargetPos=PosB, PalmTreeID = palmtrees[[i, "PalmTreeID"]]),
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
save.image(file = "~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

palmtrees %>% mutate(Start = FirstElement, End=LastElement) %>% dplyr::select(Sample, SVCaller, Chr, Start, End) %>%
write.table(file="~/Desktop/PalmTrees/Results/PalmTrees.tsv", row.names=F, col.names = T, sep = "\t", quote=F)

txptinfo %>% 
  #select(-Delly,-DellyLowQual,-Brass,-Novobreak,-Smufin,-Svaba) %>%
  write.table("~/Desktop/PalmTrees/Results/Tables/TxInPalmTrees.tsv",
              quote=F,
              sep="\t",
              row.names=F,
              col.names=T)

txptinfo %>% 
  filter(SVCaller == "Union") %>%
  distinct() %>%
  write.table("~/Desktop/PalmTrees/Results/Tables/TxInPalmTrees.tsv",
              quote=F,
              sep="\t",
              row.names=F,
              col.names=T)

