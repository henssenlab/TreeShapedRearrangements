rm(list=ls())
library(ggplot2)
library(dplyr)

load("~/Desktop/PalmTrees/Results/Figures/RelationToCSCirclesEcc/CSCircleEccOverlap_RandomSampling.Rdata")
all_all_random_regions = all_random_regions

load("~/Desktop/PalmTrees/Results/Figures/RelationToCSCirclesEc/CSCircleEcOverlap_RandomSampling.Rdata")
all_all_random_regions = rbind(all_all_random_regions, all_random_regions)

load("~/Desktop/PalmTrees/Results/Figures/RelationToWGSCirclesEcc/WGSCircleEccOverlap_RandomSampling.Rdata")
all_all_random_regions = rbind(all_all_random_regions, all_random_regions)

load("~/Desktop/PalmTrees/Results/Figures/RelationToWGSCirclesEc/WGSCircleEcOverlap_RandomSampling.Rdata")
all_all_random_regions = rbind(all_all_random_regions, all_random_regions)

#load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
#load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Circles.Rdata")

palmtree_ids_cs = palmtrees %>% filter(Sample %in% unique(cscircles_ecc$Sample)) %>% .$PalmTreeID %>% unique()
palmtree_ids_wgs = palmtrees %>% .$PalmTreeID %>% unique()

svcallers = c("Union")

for (svi in 1:length(svcallers)){
  
  thissvcaller_all_random_regions = all_all_random_regions %>% filter(SVCaller == svcallers[svi])

  thissvcaller_all_random_regions %>%
    filter(DrawID == 1, CircleCalling == "CircleSeq_ecDNA") %>%
    dplyr::select(PalmTreeID, PalmTreeOverlap, CircleCalling) %>% 
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_OverlapPalmTreeWithCSecDNA.txt"),
                quote=F, sep="\t", col.names=T, row.names=F)

  thissvcaller_all_random_regions %>%
    filter(DrawID == 1, CircleCalling == "CircleSeq_eccDNA") %>%
    dplyr::select(PalmTreeID, PalmTreeOverlap, CircleCalling) %>% 
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_OverlapPalmTreeWithCSeccDNA.txt"),
                quote=F, sep="\t", col.names=T, row.names=F)

  thissvcaller_all_random_regions %>%
    filter(DrawID == 1, CircleCalling == "WGS_ecDNA") %>%
    dplyr::select(PalmTreeID, PalmTreeOverlap, CircleCalling) %>%
    distinct() %>% 
    mutate(OverlapAtAll = PalmTreeOverlap>0.01,
           OverlapMin10Perc = PalmTreeOverlap>0.1,
           OverlapMin50Perc =  PalmTreeOverlap>0.5,
           OverlapMin75Perc =  PalmTreeOverlap>0.75,
           OverlapMin80Perc = PalmTreeOverlap>0.8,
           OverlapMin95Perc = PalmTreeOverlap>0.8) %>%
    summarise(n=n(),
              OverlapAtAll = mean(OverlapAtAll),
              OverlapMin10Perc = mean(OverlapMin10Perc),
              OverlapMin50Perc = mean(OverlapMin50Perc),
              OverlapMin75Perc = mean(OverlapMin75Perc),
              OverlapMin95Perc = mean(OverlapMin95Perc)) %>% 
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_OverlapPalmTreeWithWGSecDNA.txt"),
                quote=F, sep="\t", col.names=T, row.names=F)
    
  thissvcaller_all_random_regions %>%
      filter(DrawID == 1, CircleCalling == "WGS_eccDNA") %>%
      dplyr::select(PalmTreeID, PalmTreeOverlap, CircleCalling) %>%
      distinct() %>%
    mutate(OverlapAtAll = PalmTreeOverlap>0.01,
           OverlapMin10Perc = PalmTreeOverlap>0.1,
           OverlapMin50Perc =  PalmTreeOverlap>0.5,
           OverlapMin75Perc =  PalmTreeOverlap>0.75,
           OverlapMin80Perc = PalmTreeOverlap>0.8) %>%
      summarise(n=n(),
                OverlapAtAll = mean(OverlapAtAll),
                OverlapMin10Perc = mean(OverlapMin10Perc),
                OverlapMin50Perc = mean(OverlapMin50Perc),
                OverlapMin75Perc = mean(OverlapMin75Perc),
                OverlapMin80Perc = mean(OverlapMin80Perc)) %>% 
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_OverlapPalmTreeWithWGSeccDNA.txt"),
                quote=F, sep="\t", col.names=T, row.names=F)

MeanDistr = thissvcaller_all_random_regions %>%
    group_by(DrawID, CircleCalling) %>%
    summarise(MeanOverlap = mean(Overlap), PalmTreeOverlap = mean(PalmTreeOverlap))

  MeanDistr_cs_ecc = thissvcaller_all_random_regions %>% filter(CircleCalling=="CircleSeq_eccDNA")

  myECDF_cs_ecc = ecdf(MeanDistr %>% filter(CircleCalling=="CircleSeq_eccDNA") %>% .$MeanOverlap)
  myECDF_cs_ec = ecdf(MeanDistr %>% filter(CircleCalling=="CircleSeq_ecDNA") %>% .$MeanOverlap)
  myECDF_wgs_ecc = ecdf(MeanDistr %>% filter(CircleCalling=="WGS_eccDNA") %>% .$MeanOverlap)
  myECDF_wgs_ec = ecdf(MeanDistr %>% filter(CircleCalling=="WGS_ecDNA") %>% .$MeanOverlap)

  palmtreeoverlap_cs_ecc = MeanDistr %>% filter(DrawID == 1, CircleCalling=="CircleSeq_eccDNA") %>% .$PalmTreeOverlap %>% as.numeric()
  palmtreeoverlap_cs_ec = MeanDistr %>% filter(DrawID == 1, CircleCalling=="CircleSeq_ecDNA") %>% .$PalmTreeOverlap %>% as.numeric()
  palmtreeoverlap_wgs_ecc = MeanDistr %>% filter(DrawID == 1, CircleCalling=="WGS_eccDNA") %>% .$PalmTreeOverlap %>% as.numeric()
  palmtreeoverlap_wgs_ec = MeanDistr %>% filter(DrawID == 1, CircleCalling=="WGS_ecDNA") %>% .$PalmTreeOverlap %>% as.numeric()

  ActualPalmTreeOverlap =
    data.frame(
      "Method" = c("CircleSeq", "CircleSeq", "WGS", "WGS"),
      "CircleClass" = c("eccDNA", "ecDNA", "eccDNA", "ecDNA"),
      "MeanOverlap" = c(palmtreeoverlap_cs_ecc, palmtreeoverlap_cs_ec, palmtreeoverlap_wgs_ecc, palmtreeoverlap_wgs_ec)
    )
  
  MeanDistr %>%
    tidyr::separate(CircleCalling, into=c("Method", "CircleClass")) %>% 
  ggplot(aes(x = MeanOverlap, linetype=Method)) +
    geom_density() +
    geom_vline(data=ActualPalmTreeOverlap, aes(xintercept=MeanOverlap, linetype=Method)) +
    ggtitle(paste0("Permutation Test for Overlap\nof eccDNA and ecDNA\nwith Palm Tree Regions\n(Circle-Seq and WGS calls; ", svcallers[svi], ")")) +
    ylab("Density") +
    facet_grid(CircleClass ~.)+
    xlab(paste0("Mean Overlap with eccDNA and ecDNA \n Random Intervals (black) vs. Palm Trees (red and blue) \n", sprintf("Circle-Seq eccDNA p=%.2f\n", 1-myECDF_cs_ecc(palmtreeoverlap_cs_ecc)), sprintf("Circle-Seq ecDNA p=%.2f\n", 1-myECDF_cs_ec(palmtreeoverlap_cs_ec)), sprintf("WGS eccDNA p=%.2f\n", 1-myECDF_wgs_ecc(palmtreeoverlap_wgs_ecc)), sprintf("WGS ecDNA p=%.2f\n", 1-myECDF_wgs_ec(palmtreeoverlap_wgs_ec)))) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10),
          axis.title =  element_text(size = 10),
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) +
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToWGSandCSCircles/WGSandCSCircleOverlap_", svcallers[svi], "_PermutationTest.pdf"), height=4, width=5)

  ## Plotting
  distr = rbind(data.frame(Overlap = thissvcaller_all_random_regions$Overlap,
                           CircleCalling = thissvcaller_all_random_regions$CircleCalling,
                           isPerm = "Random Intervals"),
                data.frame(Overlap = thissvcaller_all_random_regions %>% filter(DrawID == 1) %>% .$PalmTreeOverlap,
                           CircleCalling = thissvcaller_all_random_regions %>% filter(DrawID == 1) %>% .$CircleCalling,
                           isPerm = "Palm Trees"))
  
  distr %>% 
    tidyr::separate(CircleCalling, into=c("Method", "CircleClass")) %>% 
  ggplot(aes(x=Overlap, color=isPerm, linetype=Method)) + 
    geom_density(bw=0.02, size=0.3) + 
    scale_x_continuous(limits=c(-0.2,1.2), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
    scale_color_manual(values = c("black", "red")) + 
    xlab("Relative Overlap") + 
    ylab("Density") + 
    theme_classic() +
    facet_grid(CircleClass~.) +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggtitle(paste0("Overlap of Palm Trees with \neccDNA and ecDNA \n(Circle-Seq and WGS calls; ", svcallers[svi], ")")) +
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToWGSandCSCircles/WGSandCSCircleOverlap_", svcallers[svi], "_PalmTreesVsRandomIntervals.pdf"), height=3, width=4)
}

all_all_random_regions %>% filter(DrawID == 1) %>% group_by(SVCaller, CircleCalling) %>% summarise(n=n_distinct(PalmTreeID)) %>%
  write.table(file="~/Desktop/PalmTrees/Results/Figures/RelationToWGSandCSCircles/PalmTreesPerSVCaller.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

