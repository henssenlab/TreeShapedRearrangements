rm(list=ls())
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(parallel)
cores <- detectCores()

## Load Data
setwd("~/Desktop/PalmTrees/")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Circles.Rdata")

palmtrees_onlymatched = palmtrees %>% filter(Sample %in% cscircles_ec$Sample)
txptinfo_onlymatched = txptinfo %>% filter(PalmTreeID %in% palmtrees_onlymatched$PalmTreeID)
matching_samples = unique(palmtrees_onlymatched$Sample)

# for the UCSC browser
samples = unique(cscircles_ec$Sample)
for (i in 1:length(samples)){
  cscircles_ec %>%
    filter(Sample == samples[i], CircleChr != "chrM") %>%
    dplyr::select(CircleChr, CircleStart, CircleEnd) %>%
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/forUCSCSession/PalmTreeAnalysis_CSCirclesEc_", samples[i], ".bed"),
              quote=F, sep="\t", row.names=F, col.names=F) 
}

## Compute overlap of each palm trees with ecDNA and generate matching random regions 
samples = unique(cscircles_ec$Sample)
palmtrees_onlymatched = palmtrees %>% filter(Sample %in% cscircles_ec$Sample)
palmtree_ids = unique(palmtrees_onlymatched$PalmTreeID)
txptinfo_onlymatched %>% write.table("~/Desktop/PalmTrees/Results/PalmTreeTx_OurSamples.txt", sep="\t", quote=F, row.names=F)

# for the UCSC browser
samples = unique(palmtrees_onlymatched$Sample)
for (i in 1:length(samples)){
  palmtrees_onlymatched %>%
    filter(SVCaller == "Union", Sample == samples[i]) %>%
    dplyr::select(Chr, FirstElement, LastElement) %>%
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/forUCSCSession/PalmTreeAnalysis_UnionPalmTrees_", samples[i], ".bed"),
                quote=F, sep="\t", row.names=F, col.names=F) 
}

# for the UCSC browser
samples = unique(cscircles_ec$Sample)
for (i in 1:length(samples)){
  cscircles_ec %>%
    filter(Sample == samples[i], CircleChr != "chrM") %>%
    dplyr::select(CircleChr, CircleStart, CircleEnd) %>%
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/CSCircleBeds/PalmTreeAnalysis_CSCirclesEc_", samples[i], ".bed"),
                quote=F, sep="\t", row.names=F, col.names=F) 
}

# iterate over palm trees
all_random_regions = mclapply(1:nrow(palmtrees_onlymatched), function(i){
  
  # get all circle regions within the palm tree region in this sample and compute the relative overlap for each of them
  circles_within_this_palmtree = cscircles_ec %>% 
    filter(Sample ==  palmtrees_onlymatched[[i, "Sample"]],
           CircleChr == palmtrees_onlymatched[[i, "Chr"]], 
           hasOverlap(palmtrees_onlymatched[[i, "FirstElement"]],
                      palmtrees_onlymatched[[i, "LastElement"]],
                      cscircles_ec$CircleStart,
                      cscircles_ec$CircleEnd)) %>%
    mutate(PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]],
           Sample = palmtrees_onlymatched[[i, "Sample"]], 
           SVCaller = palmtrees_onlymatched[[i, "SVCaller"]],
           FirstElement = palmtrees_onlymatched[[i, "FirstElement"]], 
           LastElement = palmtrees_onlymatched[[i, "LastElement"]],
           LengthOverlap = lengthOverlap(CircleStart, CircleEnd, palmtrees_onlymatched[[i, "FirstElement"]], palmtrees_onlymatched[[i, "LastElement"]]),
           PercentOverlap = LengthOverlap / (palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]]))
  
  # sum the relative overlaps for all merged circle regions
  if (nrow(circles_within_this_palmtree)==0){
    this_palmtree_overlap = 0
  }else{
    this_palmtree_overlap = sum(circles_within_this_palmtree$PercentOverlap)
  }
  
  # now generate 2000 random regions that match the current palm tree in size
  random_regions = sample_random_regions_broad(n=2000, segment_length = palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]])
  #random_regions = sample_random_regions_broad(n=10, segment_length = palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]])
  
  
  # compute overlap for each of the random regions
  random_regions$Overlap = 
    as.numeric(
      lapply(1:nrow(random_regions), function(j){
        circles_within_this_region = 
          cscircles_ec %>% 
          filter(Sample ==  palmtrees_onlymatched[[i, "Sample"]],
                 CircleChr == as.character(random_regions[j, "Chr"]), 
                 hasOverlap(random_regions[j, "Start"],
                            random_regions[j, "End"],
                            CircleStart,
                            CircleEnd)) %>%
          mutate(PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]],
                 Sample = palmtrees_onlymatched[[i, "Sample"]], 
                 SVCaller = palmtrees_onlymatched[[i, "SVCaller"]],
                 Start = random_regions[j, "Start"],
                 End = random_regions[j, "End"],
                 LengthOverlap = lengthOverlap(CircleStart, CircleEnd, Start, End),
                 PercentOverlap = LengthOverlap / (End-Start+1))
        if (nrow(circles_within_this_region)==0){
          this_region_overlap = 0
        }else{
          this_region_overlap = sum(circles_within_this_region$PercentOverlap)
        }
        return(this_region_overlap)
      }))
  
  # for each random region, add information about which palm tree a random region "belongs to" and the observed relative overlap this palm tree
  random_regions$DrawID = 1:nrow(random_regions)
  random_regions$PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]]
  random_regions$Sample = palmtrees_onlymatched[[i, "Sample"]]
  random_regions$SVCaller = palmtrees_onlymatched[[i, "SVCaller"]]
  random_regions$PalmTreeOverlap = this_palmtree_overlap
  
  return(random_regions)
  
},
mc.cores = cores)
all_random_regions = do.call(rbind, all_random_regions)
all_random_regions$CircleCalling = "CircleSeq_ecDNA"

save.image("~/Desktop/PalmTrees/Results/Figures/RelationToCSCirclesEc/CSCircleEcOverlap_RandomSampling.Rdata")

## For each SV caller, compare the distribution of observed vs. random overlaps
svcallers = c("Union") # Novobreak and Brass do not have >1 palmtrees for samples we have circle-seq for
for (svi in 1:length(svcallers)){
  
  # get only random regions for palm trees called by the current caller of interest
  thissvcaller_all_random_regions = all_random_regions %>% filter(SVCaller == svcallers[svi])
  
  # compute the mean relative overlap for all palm tree regions and for each of the 500 sets of matched random regions
  MeanDistr = thissvcaller_all_random_regions %>%
    group_by(DrawID) %>%
    summarise(MeanOverlap = mean(Overlap), PalmTreeOverlap = mean(PalmTreeOverlap))
  
  # plot the distribution of mean overlap for the 500 sets of matched random regions and compare to the mean overlap of the set of actual palm trees
  myECDF = ecdf(MeanDistr$MeanOverlap)
  ggplot(data=MeanDistr, aes(x = MeanOverlap)) + 
    geom_density() + 
    geom_vline(xintercept = MeanDistr[[1,"PalmTreeOverlap"]], color="red") + 
    ggtitle(paste0("Test for ecDNA within Palm Tree Regions (CS calls; ", svcallers[svi], ")")) +
    ylab("Density") +
    xlab(paste0("Mean Overlap with ecDNA \n Random Intervals (black) vs. Palm Trees (red) \n", sprintf("p=%.2f", 1-myECDF(MeanDistr[[1,"PalmTreeOverlap"]])))) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToCSCirclesEc/CSCircleEcOverlap_", svcallers[svi], "_PermutationTest.pdf"), height=4, width=6)
  
  
  # plot the distribution of relative overlaps in the acutal set of palm trees vs. overlaps for all randomly generated matched regions
  distr = rbind(data.frame(Overlap = thissvcaller_all_random_regions$Overlap,
                           isPerm = "Random Intervals"),
                data.frame(Overlap = thissvcaller_all_random_regions %>% filter(DrawID == 1) %>% .$PalmTreeOverlap,
                           isPerm = "Palm Trees"))
  
  ggplot(data=distr, aes(x=Overlap, color=isPerm)) + 
    geom_density(bw=0.02, size=0.3) + 
    scale_x_continuous(limits=c(-0.2,1.2), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
    scale_color_manual(values = c("black", "red")) + 
    xlab("Overlap with ecDNA") + 
    ylab("Density") + 
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggtitle(paste0("Overlap of Palm Trees with ecDNA (CS calls; ", svcallers[svi], ")")) +
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToCSCirclesEc/CSCircleEcOverlap_", svcallers[svi], "_PalmTreesVsRandomIntervals.pdf"), height=4, width=6)
}


