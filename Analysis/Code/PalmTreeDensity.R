rm(list = ls())
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
source('~/Desktop/PalmTrees/Analysis/Code/ggGrandLinear.R')

# The number of all samples that went through SV calling 
# (irrespective of whether the sample was excluded because of quality issues
# or whether any SV was called)
nsamples = 97

svcallers = unique(palmtrees$SVCaller)
hg19 <- BSgenome.Hsapiens.UCSC.hg19
bins = tileGenome(seqinfo(hg19), tilewidth=10000, cut.last.tile.in.chrom=TRUE)
seqlevels(bins, pruning.mode="coarse") = seqlevels(bins)[1:23]

for (i in 1:length(svcallers)){
  
  PalmTreeCoverage = 
    coverage(palmtrees_gr[palmtrees_gr$SVCaller == svcallers[i]]) %>%
    binnedAverage(bins, ., "Coverage") %>% 
    data.frame() %>%
    as_tibble()
  
  PalmTreeCoverage$Recurrence = ceiling(PalmTreeCoverage$Coverage) / 97
  
  # Plot Histogram 
  f=ggGrandLinear(PalmTreeCoverage$seqnames, PalmTreeCoverage$start, 100*PalmTreeCoverage$Recurrence) + 
    geom_line() +
    ylab("Percent Samples") +
    ggtitle(paste0("Palm Tree Recurrence (", svcallers[i], ")"))
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/", svcallers[i], "PalmTreeRecurence_10kbBins.eps"), plot=f, width=9, height=3)

}
