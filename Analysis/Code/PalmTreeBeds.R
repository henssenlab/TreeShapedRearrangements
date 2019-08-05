library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(biomaRt)

setwd("~/Desktop/PalmTrees/")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

bed = bind_cols(palmtrees[,"Chr"], palmtrees[,"FirstElement"], palmtrees[,"LastElement"])
samples = unique(as.character(palmtrees$Sample))
svcallers = unique(as.character(palmtrees$SVCaller))
for (sai in 1:length(samples)){
  for (svi in 1:length(svcallers)){
    this_sample_this_svcaller = palmtrees %>% filter(Sample == samples[sai], SVCaller == svcallers[svi])
    if (nrow(this_sample_this_svcaller)>0){
    this_bed = bind_cols(this_sample_this_svcaller[,"Chr"], this_sample_this_svcaller[,"FirstElement"], this_sample_this_svcaller[,"LastElement"])
    write.table(this_bed, 
                file=paste0("~/Desktop/PalmTrees/Results/PalmTreeBeds/", samples[sai], "_",svcallers[svi], ".bed"),
                col.names=FALSE,
                row.names=FALSE, 
                sep="\t",
                quote=FALSE)
    }
  }
}

# mycols = rand_color(n=100, luminosity = "dark")
# circos.clear()
# circos.par("start.degree" = 90)
# #circos.initializeWithIdeogram(species = "hg19")
# circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
# text(0, 0, "", cex = 1) 
# palmtrees_bed = bed
# diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
# if (sum(diff<5000000)>0){
#   palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff)/2
#   palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff)/2
# }
# circos.genomicTrack(as.data.frame(palmtrees_bed),
#                     panel.fun = function(region, value, ...){
#                       circos.genomicRect(region, value, col="steelblue", border=NA, ...)
#                     },
#                     ylim=c(0,0.2), track.height=0.33*circos.par("track.height"))
# 




