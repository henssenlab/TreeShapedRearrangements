rm(list=ls())
library(dplyr)
library(ggplot2)
setwd("~/Desktop/PalmTrees/")
source("~/Desktop/PalmTrees/Analysis/Code/ParseCopyNumberData.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")

# this actually does not filter out anything
palmtrees = palmtrees %>% filter(Sample %in% ascat_cnv$Sample)
palmtrees$MeanTumorCN = NA
palmtrees$MaxTumorCN = NA
palmtrees$MeanNormalCN = NA

for (i in 1:nrow(palmtrees)){
  which_sample = palmtrees[[i, "Sample"]]
  chr = palmtrees[[i, "Chr"]]
  start = palmtrees[[i, "FirstElement"]]
  end = palmtrees[[i, "LastElement"]]
  palmtree_length = end-start+1
  
  this_cnv_data = ascat_cnv %>% filter(Sample == which_sample,
                                       Chr == chr,
                                       hasOverlap(start, end, ascat_cnv$Start, ascat_cnv$End))
  if (nrow(this_cnv_data) == 0){
    this_cnv_data = data.frame(NormalTotalCopyNumber = 2, TumorTotalCopyNumber = 2, TumorFC = 1, OverlapLength = palmtree_length, OverlapLengthPercent = 1)
  }else{
    this_cnv_data$OverlapLength = lengthOverlap(start, end, this_cnv_data$Start, this_cnv_data$End)
    this_cnv_data$OverlapLengthPercent = lengthOverlap(start, end, this_cnv_data$Start, this_cnv_data$End) / palmtree_length
  }
  
  palmtrees[[i, "MeanTumorCN"]] = sum(this_cnv_data$TumorTotalCopyNumber * this_cnv_data$OverlapLengthPercent) + 2 * (1-sum(this_cnv_data$OverlapLengthPercent))
  palmtrees[[i, "MeanNormalCN"]] = sum(this_cnv_data$NormalTotalCopyNumber * this_cnv_data$OverlapLengthPercent) + 2 * (1-sum(this_cnv_data$OverlapLengthPercent))
  palmtrees[[i, "MaxTumorCN"]] = max(this_cnv_data$TumorTotalCopyNumber)
  
}

palmtrees %>%
  mutate(isAmplified = MeanTumorCN>=9,
         isGained = MeanTumorCN<9 & MeanTumorCN>=3,
         isLost = MeanTumorCN<=1,
         isCopyNumberNeutral = MeanTumorCN>1 & MeanTumorCN<3) %>% 
  group_by(SVCaller) %>% 
  summarise(
    NumberOfPalmTrees = n(),
    nAmplified = sum(isAmplified),
    PercentAmplified = mean(isAmplified),
    nGained = sum(isGained),
    PercentGained = mean(isGained),
    nCNNeutral = sum(isCopyNumberNeutral),
    PercentCNNeutral = mean(isCopyNumberNeutral),
    nLost = sum(isLost),
    PercentLost = mean(isLost)
  ) %>%
  write.table("~/Desktop/PalmTrees/Results/Tables/PalmTreesVsCopyNumber.txt", quote=F, col.names=T, row.names=F, sep="\t")



# m=lm(ElementsPerMb ~ MeanTumorCN, data=palmtrees %>% filter(ElementsPerMb<100, SVCaller == "Union"))
# palmtrees %>%
#   filter(ElementsPerMb<100, SVCaller == "Union") %>%
#   ggplot(aes(x=MeanTumorCN, y=ElementsPerMb)) +
#   geom_smooth(method="lm", se=F, color="steelblue") +
#   geom_point() +
#   theme_kons1() +
#   xlab("Palm Tree Interval Copy Number") +
#   ylab("Rearrangements per Mb") +
#   ggtitle("Union Palm Trees") +
#   annotate(geom="text", x=Inf, y=-Inf, label=sprintf("p=%.e", summary(m)$coefficients[,4][2]), hjust=1, vjust=-0.5) +
#   ggsave("~/Desktop/PalmTrees/Results/Figures/RelationToCopyNumber/Union_PalmTreeTxVsCopyNumber.pdf", height=3, width=3)
#
# m=lm(ElementsPerMb ~ MeanTumorCN, data=palmtrees %>% filter(ElementsPerMb<100, SVCaller == "Delly"))
# palmtrees %>%
#   filter(ElementsPerMb<100, SVCaller == "Delly") %>%
#   ggplot(aes(x=MeanTumorCN, y=ElementsPerMb)) +
#   geom_smooth(method="lm", se=F, color="steelblue") +
#   geom_point() +
#   theme_kons1() +
#   xlab("Palm Tree Interval Copy Number") +
#   ylab("Rearrangements per Mb") +
#   ggtitle("Delly Palm Trees") +
#   annotate(geom="text", x=Inf, y=-Inf, label=sprintf("p=%.e", summary(m)$coefficients[,4][2]), hjust=1, vjust=-0.5) +
#   ggsave("~/Desktop/PalmTrees/Results/Figures/RelationToCopyNumber/Delly_PalmTreeTxVsCopyNumber.pdf", height=3, width=3)
#
# m=lm(ElementsPerMb ~ MeanTumorCN, data=palmtrees %>% filter(ElementsPerMb<100, SVCaller == "Smufin"))
# palmtrees %>%
#   filter(ElementsPerMb<100, SVCaller == "Smufin") %>%
#   ggplot(aes(x=MeanTumorCN, y=ElementsPerMb)) +
#   geom_smooth(method="lm", se=F, color="steelblue") +
#   geom_point() +
#   theme_kons1() +
#   xlab("Palm Tree Interval Copy Number") +
#   ylab("Rearrangements per Mb") +
#   ggtitle("Smufin Palm Trees") +
#   annotate(geom="text", x=Inf, y=-Inf, label=sprintf("p=%.e", summary(m)$coefficients[,4][2]), hjust=1, vjust=-0.5) +
#   ggsave("~/Desktop/PalmTrees/Results/Figures/RelationToCopyNumber/Smufin_PalmTreeTxVsCopyNumber.pdf", height=3, width=3)
#
# palmtrees %>%
#   group_by(SVCaller) %>%
#   summarise(MeanCopyNumber = mean(MeanTumorCN))

