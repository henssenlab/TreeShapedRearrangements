rm(list=ls())
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(parallel)
cores <- detectCores()

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")

system("rm  ~/Desktop/PalmTrees/Data/WGSCircles/*.bed; cp ~/Desktop/eccDNADataAnalysis/*Nov2018/*_tumor_*Sort.bed ~/Desktop/PalmTrees/Data/WGSCircles/")
system("for i in ~/Desktop/PalmTrees/Data/WGSCircles/*Sort.bed; do awk '{print $0 \"\\t\" FILENAME}' $i > ${i}.tmp ; mv ${i}.tmp ${i} ; done")
system("awk FNR-1 ~/Desktop/PalmTrees/Data/WGSCircles/*Sort.bed > ~/Desktop/PalmTrees/Data/WGSCircles/AllSamples_Sort.bed")

wgscircles = read.table("~/Desktop/PalmTrees/Data/WGSCircles/AllSamples_Sort.bed", sep="\t", header=F)
colnames(wgscircles) = c("CircleChr", "CircleStart", "CircleEnd", "DontKnow1", "CircleLength", "DontKnow2", "DontKnow3", "DontKnow4", "DontKnow5", "Sample")
wgscircles$CircleChr = as.character(paste0("chr", as.character(wgscircles$CircleChr)))
wgscircles$Sample = lapply(wgscircles$Sample, function (s) gsub("_.*", "", gsub(".*WGSCircles/", "", s)))

# Merge those intervals
merged_wgscircles = list()
samples = unique(wgscircles$Sample)

merged_wgscircles = lapply(samples, function (sample){
  this_circles = wgscircles %>% filter(Sample == sample)
  if (nrow(this_circles)==0) return(data.frame(MergedCirclesChr = NA, MergedCirclesStart = NA, MergedCirclesEnd = NA, MergedCirclesLength = NA, Sample = sample))
  this_circles = makeGRangesFromDataFrame(df = this_circles,
                                          seqnames.field = "CircleChr",
                                          start.field = "CircleStart",
                                          end.field = "CircleEnd",
                                          keep.extra.columns = T)
  this_circles = data.frame(reduce(this_circles)) %>%
    dplyr::select(seqnames, start, end, width)
  colnames(this_circles) = c("MergedCirclesChr", "MergedCirclesStart", "MergedCirclesEnd", "MergedCirclesLength")
  this_circles$Sample = sample
  return(this_circles)
})
merged_wgscircles = do.call(rbind, merged_wgscircles)
merged_wgscircles = merged_wgscircles %>% filter(isdefchrom(MergedCirclesChr)) %>% droplevels()

## Intersect with palmtree
palmtrees_onlymatched = palmtrees %>% filter(Sample %in% merged_wgscircles$Sample)
palmtree_ids = unique(palmtrees_onlymatched$PalmTreeID)

all_random_regions = mclapply(1:nrow(palmtrees_onlymatched), function(i){
  
  circles_within_this_palmtree = merged_wgscircles %>% 
    filter(Sample ==  palmtrees_onlymatched[[i, "Sample"]],
           MergedCirclesChr == palmtrees_onlymatched[[i, "Chr"]], 
           hasOverlap(palmtrees_onlymatched[[i, "FirstElement"]],
                      palmtrees_onlymatched[[i, "LastElement"]],
                      merged_wgscircles$MergedCirclesStart,
                      merged_wgscircles$MergedCirclesEnd)) %>%
    mutate(PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]],
           Sample = palmtrees_onlymatched[[i, "Sample"]], 
           SVCaller = palmtrees_onlymatched[[i, "SVCaller"]],
           FirstElement = palmtrees_onlymatched[[i, "FirstElement"]], 
           LastElement = palmtrees_onlymatched[[i, "LastElement"]],
           LengthOverlap = lengthOverlap(MergedCirclesStart, MergedCirclesEnd, palmtrees_onlymatched[[i, "FirstElement"]], palmtrees_onlymatched[[i, "LastElement"]]),
           PercentOverlap = LengthOverlap / (palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]]))
  
  if (nrow(circles_within_this_palmtree)==0){
    this_palmtree_overlap = 0
  }else{
    this_palmtree_overlap = sum(circles_within_this_palmtree$PercentOverlap)
  }
  
  # now generate a bunch of random regions
  random_regions = sample_random_regions(n=500, segment_length = palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]])
  
  random_regions$Overlap = 
    as.numeric(
      lapply(1:nrow(random_regions), function(j){
        circles_within_this_region = 
          merged_wgscircles %>% 
          filter(Sample ==  palmtrees_onlymatched[[i, "Sample"]],
                 MergedCirclesChr == as.character(random_regions[j, "Chr"]), 
                 hasOverlap(random_regions[j, "Start"],
                            random_regions[j, "End"],
                            MergedCirclesStart,
                            MergedCirclesEnd)) %>%
          mutate(PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]],
                 Sample = palmtrees_onlymatched[[i, "Sample"]], 
                 SVCaller = palmtrees_onlymatched[[i, "SVCaller"]],
                 Start = random_regions[j, "Start"],
                 End = random_regions[j, "End"],
                 LengthOverlap = lengthOverlap(MergedCirclesStart, MergedCirclesEnd, Start, End),
                 PercentOverlap = LengthOverlap / (End-Start+1))
        if (nrow(circles_within_this_region)==0){
          this_region_overlap = 0
        }else{
          this_region_overlap = sum(circles_within_this_region$PercentOverlap)
        }
        return(this_region_overlap)
      }))
  
  random_regions$DrawID = 1:nrow(random_regions)
  random_regions$PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]]
  random_regions$Sample = palmtrees_onlymatched[[i, "Sample"]]
  random_regions$SVCaller = palmtrees_onlymatched[[i, "SVCaller"]]
  random_regions$PalmTreeOverlap = this_palmtree_overlap
  
  return(random_regions)
  
},
mc.cores = cores)

all_random_regions_wgs = do.call(rbind, all_random_regions)

save.image("~/Desktop/PalmTrees/Results/Figures/RelationToWGSCircles/WGSCircleOverlap_RandomSampling500_181213.Rdata")

## Intersect with palmtree
match_names = read.table("~/Desktop/PalmTrees/Data/MatchNames.csv", sep=";", header=T)
match_names$cologne=as.character(match_names$matchID.Cologne.)
match_names$name=as.character(match_names$WGS_ID)
meta = read.table("~/Desktop/eccDNADataAnalysis/Richard\ tables/SamplesKonstantinDataAnalysis/metaeccdna.tsv", header=T)
meta$name = as.character(meta$name)
make_fname_thr3 = function (sample){
  if (substr(sample,1,3)=="BER"){
    paste0("~/Desktop/eccDNADataAnalysis/Richard\ tables/SamplesKonstantinDataAnalysis/Homer_findPeaks_histMode_", sample, ".hg19.sorted.RmDup.regions_merge3kb_SAF_CrcReadCounts_Thresh3.txt")
  }else{
    paste0("~/Desktop/eccDNADataAnalysis/Richard\ tables/SamplesKonstantinDataAnalysis/Homer_findPeaks_histMode_", sample, ".hg19.sorted.RmDup.regions_merge3kb_toSAF_", sample, ".hg19.sorted.RmDup.CircPairTestF49.QFilt20.merged_ftCounts_Thresh3.txt")
  }
}
# filter palmtree data for samples we also have eccDNA data for
matches = inner_join(meta, match_names, by="name") %>% dplyr::select(fname_infix, ecc_ID, name, cologne)
palmtrees_onlymatched = palmtrees %>% filter(Sample %in% matches$cologne)

all_random_regions = mclapply(1:nrow(palmtrees_onlymatched), function(i){
  
  eccDNA = read.table(matches %>% filter(cologne == palmtrees_onlymatched[[i,"Sample"]]) %>% .$fname_infix %>% as.character %>% make_fname_thr3())
  
  circles_within_this_palmtree = merged_wgscircles %>% 
    filter(Sample ==  palmtrees_onlymatched[[i, "Sample"]],
           MergedCirclesChr == palmtrees_onlymatched[[i, "Chr"]], 
           hasOverlap(palmtrees_onlymatched[[i, "FirstElement"]],
                      palmtrees_onlymatched[[i, "LastElement"]],
                      merged_wgscircles$MergedCirclesStart,
                      merged_wgscircles$MergedCirclesEnd)) %>%
    mutate(PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]],
           Sample = palmtrees_onlymatched[[i, "Sample"]], 
           SVCaller = palmtrees_onlymatched[[i, "SVCaller"]],
           FirstElement = palmtrees_onlymatched[[i, "FirstElement"]], 
           LastElement = palmtrees_onlymatched[[i, "LastElement"]],
           LengthOverlap = lengthOverlap(MergedCirclesStart, MergedCirclesEnd, palmtrees_onlymatched[[i, "FirstElement"]], palmtrees_onlymatched[[i, "LastElement"]]),
           PercentOverlap = LengthOverlap / (palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]]))
  
  if (nrow(circles_within_this_palmtree)==0){
    this_palmtree_overlap = 0
  }else{
    this_palmtree_overlap = sum(circles_within_this_palmtree$PercentOverlap)
  }
  
  # now generate a bunch of random regions
  random_regions = sample_random_regions(n=500, segment_length = palmtrees_onlymatched[[i, "LastElement"]]-palmtrees_onlymatched[[i, "FirstElement"]])
  
  random_regions$Overlap = 
    as.numeric(
      lapply(1:nrow(random_regions), function(j){
        circles_within_this_region = 
          merged_wgscircles %>% 
          filter(Sample ==  palmtrees_onlymatched[[i, "Sample"]],
                 MergedCirclesChr == as.character(random_regions[j, "Chr"]), 
                 hasOverlap(random_regions[j, "Start"],
                            random_regions[j, "End"],
                            MergedCirclesStart,
                            MergedCirclesEnd)) %>%
          mutate(PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]],
                 Sample = palmtrees_onlymatched[[i, "Sample"]], 
                 SVCaller = palmtrees_onlymatched[[i, "SVCaller"]],
                 Start = random_regions[j, "Start"],
                 End = random_regions[j, "End"],
                 LengthOverlap = lengthOverlap(MergedCirclesStart, MergedCirclesEnd, Start, End),
                 PercentOverlap = LengthOverlap / (End-Start+1))
        if (nrow(circles_within_this_region)==0){
          this_region_overlap = 0
        }else{
          this_region_overlap = sum(circles_within_this_region$PercentOverlap)
        }
        return(this_region_overlap)
      }))
  
  random_regions$DrawID = 1:nrow(random_regions)
  random_regions$PalmTreeID = palmtrees_onlymatched[[i, "PalmTreeID"]]
  random_regions$Sample = palmtrees_onlymatched[[i, "Sample"]]
  random_regions$SVCaller = palmtrees_onlymatched[[i, "SVCaller"]]
  random_regions$PalmTreeOverlap = this_palmtree_overlap
  
  return(random_regions)
  
},
mc.cores = cores)

all_random_regions_ecc = do.call(rbind, all_random_regions)

####

svcallers = c("Smufin", "Delly", "Svaba", "Brass", "AtLeastTwo", "Novobreak")

for (svi in 1:length(svcallers)){
  
  thissvcaller_all_random_regions = all_random_regions %>% filter(SVCaller == svcallers[svi])
  
  MeanDistr = thissvcaller_all_random_regions %>%
    group_by(DrawID) %>%
    summarise(MeanOverlap = mean(Overlap), PalmTreeOverlap = mean(PalmTreeOverlap))
  
  myECDF = ecdf(MeanDistr$MeanOverlap)
  
  ggplot(data=MeanDistr, aes(x = MeanOverlap)) + 
    geom_density() + 
    geom_vline(xintercept = MeanDistr[[1,"PalmTreeOverlap"]], color="red") + 
    ggtitle(paste0("Test for eccDNA within Palm Tree Regions (WGS calls; ", svcallers[svi], ")")) +
    ylab("Density") +
    xlab(paste0("Mean Overlap with eccDNA \n Random Intervals (black) vs. Palm Trees (red) \n", sprintf("p=%.2f", 1-myECDF(MeanDistr[[1,"PalmTreeOverlap"]])))) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToWGSCircles/WGSCircleOverlap_", svcallers[svi], "_PermutationTest.pdf"), height=4, width=6)
  
  
  ## Plotting
  distr = rbind(data.frame(Overlap = thissvcaller_all_random_regions$Overlap,
                           isPerm = "Random Intervals"),
                data.frame(Overlap = thissvcaller_all_random_regions %>% filter(DrawID == 1) %>% .$PalmTreeOverlap,
                           isPerm = "Palm Trees"))
  
  ggplot(data=distr, aes(x=Overlap, color=isPerm)) + 
    geom_density(bw=0.02, size=0.3) + 
    scale_x_continuous(limits=c(-0.2,1.2), breaks=c(0,0.2,0.4,0.6,0.8,1)) +
    scale_color_manual(values = c("black", "red")) + 
    xlab("Overlap with eccDNA") + 
    ylab("Density") + 
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggtitle(paste0("Overlap of Palm Trees with eccDNA (WGS calls; ", svcallers[svi], ")")) +
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToWGSCircles/WGSCircleOverlap_", svcallers[svi], "_PalmTreesVsRandomIntervals.pdf"), height=4, width=6)
}

