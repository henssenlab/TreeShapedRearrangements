rm(list = ls())

library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
source('~/Desktop/PalmTrees/Analysis/Code/ggGrandLinear.R')

#################
callinginfo = read.csv("~/Desktop/PalmTrees/Data/EliasMail181205/SampleCalls.csv", header=T, sep=";")
callinginfo$AtLeastTwo = 1 # that is not a good idea, but maybe the best
callinginfo$Delly = 1 # that is not a good idea, but maybe the best
callinginfo$Smufin = 1 # that is not a good idea, but maybe the best
callinginfo$Novobreak = 1 # that is not a good idea, but maybe the best
callinginfo$Brass = 1 # that is not a good idea, but maybe the best
callinginfo$Union = 1 # that is not a good idea, but maybe the best

#################

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

PalmTreeCoverage$Recurrence = ceiling(PalmTreeCoverage$Coverage) / sum(callinginfo[,svcallers[i]])

f=ggGrandLinear(PalmTreeCoverage$seqnames, PalmTreeCoverage$start, 100*PalmTreeCoverage$Recurrence) + 
  geom_smooth(method="gam", formula=y ~ s(x, bs="cs"), se=FALSE, method.args=list(family=poisson)) +
  ylab("Percent Samples") +
  ggtitle(paste0("Palm Tree Recurrence (", svcallers[i], ")"))
ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/", svcallers[i], "PalmTreeRecurrence_10kbBins_Smooth.eps"), plot=f, width=9, height=3)

f=ggGrandLinear(PalmTreeCoverage$seqnames, PalmTreeCoverage$start, 100*PalmTreeCoverage$Recurrence) + 
  geom_line() +
  ylab("Percent Samples") +
  ggtitle(paste0("Palm Tree Recurrence (", svcallers[i], ")"))
ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/", svcallers[i], "PalmTreeRecurence_10kbBins.eps"), plot=f, width=9, height=3)

}

# #################
# #seqlevels(palmtrees_gr, pruning.mode="coarse") = seqlevels(bins)
# 
# PalmTreeCoverage = coverage(palmtrees_gr) %>%
#   binnedAverage(bins, ., "Coverage") %>% 
#   data.frame() %>%
#   as_tibble()
# 
# 
# ggGrandLinear(PalmTreeCoverage$seqnames, PalmTreeCoverage$start, PalmTreeCoverage$Coverage) + 
#   geom_point()
# ggsave("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/PalmTreeCoverageTileWidth10kbManhattan.pdf")
# 
# ggGrandLinear(PalmTreeCoverage$seqnames, PalmTreeCoverage$start, PalmTreeCoverage$Coverage) + 
#   geom_smooth(se=F, span=0.4)
# ggsave("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/PalmTreeCoverageTileWidth10kbDensity.pdf")
# 
# ggGrandLinearGenomeWideDensity(palmtrees$Chr, palmtrees$FirstElement, bw=1)
# ggsave("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/PalmTreeSingleElementsDensity.pdf")
# 
# ggGrandLinearGenomeWideHistogram(palmtrees, "Chr", "FirstElement", binwidthMb=15) 
# ggsave("~/Desktop/PalmTrees/Results/Figures/PalmTreeDistribution/PalmTreeSingleElementsHistogramBinwidth15Mb.pdf")
# 
