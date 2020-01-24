rm(list=ls())

library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")

source("~/Desktop/PalmTrees/Analysis/Code/ParseAllSV.R")

datapath = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/"
samples = dir(datapath)
samples = samples[!grepl("excluded", samples)]

svaba = list()
for (i in 1:length(samples)){
  svaba[[i]] = parse_svaba_allsv(paste0(datapath, samples[i], "/svaba/svaba-indels.PASS.vcf"), paste0(datapath, samples[i],"/svaba/svaba-sv-fromindels.PASS.vcf"), paste0(datapath, samples[i],"/svaba/svaba-sv.PASS.vcf"), samples[i])
}
svaba = do.call(rbind, svaba) 
svaba = svaba %>% filter(Class == "sv", Filter=="PASS") %>% as_tibble()
svaba$Sample = gsub("CB", "NB", svaba$Sample)
svaba = svaba %>% filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
svaba_allsv = svaba

svaba_ins_hom = svaba_allsv %>% 
  select(Sample, BNDPairID, Insertion, InsertionLength, Homology, HomologyLength,
         ChrA, PosA, ChrB, PosB, DirectionA, DirectionB) %>%
  mutate(SampleBNDPairID = paste0(Sample, "_", BNDPairID))

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/SvabaSV.Rdata")

