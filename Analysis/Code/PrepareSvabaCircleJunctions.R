rm(list=ls())

library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")
source("~/Desktop/PalmTrees/Analysis/Code/ParseSvabaCircleJunctionVCF.R")

matchnames = read.table("~/Desktop/PalmTrees/Data/MatchNames.csv", header=T, sep=";")
colnames(matchnames) = c("eccID", "Sample", "WGS_ID")

sv_files = list.files(path="~/Desktop/PalmTrees/Data/CircleSeq_CircleReads_Svaba/", pattern="*.bam.svaba.sv.vcf")
sample_names = c("Kelly", "SKNSH", "VH7", "SKNFI", "IMR575", "SHSY5Y", "SKNAS", "GIMEN", "NB69",
                 "CHP212", "RPE", "RPE-GFP", "RPE-GFP-PGBD5", "A68277", "B3359564", "18149",
                 "21533", "25225", "20NBLPt", "21NBLPt", "RH4", "SHEP", "NGP", "MA6", "VM1", 
                 "26891", "23269", "23484", "23570", "24912", "18800", "UW228", "20285", "20420", "20494", "20477", "21258", "21022", 
                 "TE381T", "G401", "RD", "RH30", "T174")
info = data.frame(
  "fname" = sv_files,
  "WGS_ID" = sample_names
)

# translate WGS_ID to actual sample names
info = info %>% left_join(matchnames)

# filter all non-patients out, i.e. those who do not have a sample name in matchnames
info = info %>% filter(!is.na(Sample))

svaba_junctions=list()
for (i in 1:nrow(info)){
  print(as.character(info[i, "Sample"]))
  svaba_junctions[[i]] = parse_svaba_circlejunction_vcf(paste0("~/Desktop/PalmTrees/Data/CircleSeq_CircleReads_Svaba/",info[i, "fname"]), as.character(info[i, "Sample"]))
}
svaba_junctions = do.call(rbind, svaba_junctions)

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/SvabaCircleJunctions.Rdata")

