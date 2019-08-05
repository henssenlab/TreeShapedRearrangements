rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/GeneCopyNumber.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/PalmTreeXExpressionData.Rdata")

gene_cn = gene_cn %>% mutate(ASCAT_TumorTotalCopyNumber = CopyNumber) %>% dplyr::select(Sample, Gene, ASCAT_TumorTotalCopyNumber)

cosmic = read.table("~/Desktop/PalmTrees/Data/COSMIC_Census_allTue Dec 11 07_52_32 2018.tsv", header=T, sep="\t")

completedata %>%
  dplyr::select(Sample, SVCaller,Gene, PalmTreeID, CoordDistFromBreakpoint, TPM, FoldChange, ModifiedZScore) %>%
  full_join(data.frame(Gene=cosmic$Gene.Symbol, COSMIC_Annotation=cosmic$Role.in.Cancer)) %>%
  filter(SVCaller == "Union", (CoordDistFromBreakpoint<=500000 & CoordDistFromBreakpoint>=-500000), (ModifiedZScore >= 2 | ModifiedZScore < -1)) %>% 
  left_join(gene_cn, by=c("Sample", "Gene")) %>%
  write.table("~/Desktop/PalmTrees/Results/Tables/Union_DeregulatedGenesAroundTargetSites.txt", sep="\t",quote=F, col.names=T,row.names=F)
