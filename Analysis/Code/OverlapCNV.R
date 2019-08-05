library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Desktop/PalmTrees/")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")

# this actually does not filter out anything
palmtrees = palmtrees %>% filter(Sample %in% ascat_cnv$Sample)
palmtrees$MeanTumorCN = NA
palmtrees$MeanNormalCN = NA

all_random_regions_tumor_mean_cn = list()
all_random_regions_normal_mean_cn = list()
all_random_regions_which_palmtree = list()
all_random_regions_svcaller = list()
all_random_regions_sample = list()
all_random_regions_draw_id = list()

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
  palmtrees[[i, "MeanTumorCN"]] = mean(this_cnv_data$TumorTotalCopyNumber) * sum(this_cnv_data$OverlapLengthPercent) + 2 * (1-sum(this_cnv_data$OverlapLengthPercent))
  palmtrees[[i, "MeanNormalCN"]] = mean(this_cnv_data$NormalTotalCopyNumber) * sum(this_cnv_data$OverlapLengthPercent) + 2 * (1-sum(this_cnv_data$OverlapLengthPercent))

  random_regions = sample_random_regions(n=500, segment_length = palmtree_length)
  random_regions_tumor_mean_cn = NA*(1:nrow(random_regions))
  random_regions_normal_mean_cn = NA*(1:nrow(random_regions))
  random_regions_which_palmtree = NA*(1:nrow(random_regions))
  random_regions_svcaller = NA*(1:nrow(random_regions))
  random_regions_sample = NA*(1:nrow(random_regions))
  random_regions_draw_id = NA*(1:nrow(random_regions))
  for (j in 1:nrow(random_regions)){
    this_random_region_cnv_data = ascat_cnv %>% filter(Sample == which_sample,
                                                       Chr == chr,
                                                       hasOverlap(random_regions[j,"Start"], random_regions[j,"End"], ascat_cnv$Start, ascat_cnv$End))
    if (nrow(this_random_region_cnv_data) == 0){
      this_random_region_cnv_data = data.frame(NormalTotalCopyNumber = 2, TumorTotalCopyNumber = 2, OverlapLength = palmtree_length, OverlapLengthPercent = 1)
    }else{
      this_random_region_cnv_data$OverlapLength = lengthOverlap(random_regions[j,"Start"], random_regions[j,"End"], this_random_region_cnv_data$Start, this_random_region_cnv_data$End)
      this_random_region_cnv_data$OverlapLengthPercent = lengthOverlap(random_regions[j,"Start"], random_regions[j,"End"], this_random_region_cnv_data$Start, this_random_region_cnv_data$End) / palmtree_length
    }
    random_regions_tumor_mean_cn[j] = mean(this_random_region_cnv_data$TumorTotalCopyNumber) * sum(this_random_region_cnv_data$OverlapLengthPercent) + 2 * (1-sum(this_random_region_cnv_data$OverlapLengthPercent))
    random_regions_normal_mean_cn[j] = mean(this_random_region_cnv_data$NormalTotalCopyNumber) * sum(this_random_region_cnv_data$OverlapLengthPercent) + 2 * (1-sum(this_random_region_cnv_data$OverlapLengthPercent))
    random_regions_which_palmtree[j] = palmtrees[[i, "PalmTreeID"]]
    random_regions_svcaller[j] = palmtrees[[i, "SVCaller"]]
    random_regions_sample[j] = palmtrees[[i, "Sample"]]
    random_regions_draw_id[j] = j
    
  }
  
  all_random_regions_tumor_mean_cn[[i]] = random_regions_tumor_mean_cn
  all_random_regions_normal_mean_cn[[i]] = random_regions_normal_mean_cn
  all_random_regions_which_palmtree[[i]] = random_regions_which_palmtree
  all_random_regions_svcaller[[i]] = random_regions_svcaller
  all_random_regions_sample[[i]] = random_regions_sample
  all_random_regions_draw_id[[i]] = random_regions_draw_id
  
  }



## Perm test
RandomDraws = data.frame(TumorCN = unlist(all_random_regions_tumor_mean_cn),
                         NormalCN = unlist(all_random_regions_normal_mean_cn),
                         PalmTreeID = unlist(all_random_regions_which_palmtree),
                         SVCaller = unlist(all_random_regions_svcaller),
                         Sample = unlist(all_random_regions_sample),
                         DrawIndex = unlist(all_random_regions_draw_id))
RandomDraws$logFC = log2( (RandomDraws$TumorCN+1) / (RandomDraws$NormalCN+1) )
palmtrees$logFC = log2( (palmtrees$MeanTumorCN+1) / (palmtrees$MeanNormalCN+1) )

save.image("~/Desktop/PalmTrees/Results/Figures/RelationToCopyNumber/RandomSampling500_181213.Rdata")

svcallers = c("Smufin", "Delly", "Brass", "Svaba", "AtLeastTwo", "Novobreak")

for (svi in 1:length(svcallers)){
  
  this_svcaller_RandomDraws = RandomDraws %>% filter(SVCaller == svcallers[svi])
  this_svcaller_palmtrees = palmtrees %>% filter(SVCaller == svcallers[svi])
  
  MeanDistr = this_svcaller_RandomDraws %>%
    group_by(DrawIndex) %>%
    summarise(MeanLogFC = mean(logFC)) %>% 
    arrange(desc(MeanLogFC))
  
  myECDF = ecdf(MeanDistr$MeanLogFC)
  
  ggplot(data=MeanDistr, aes(x = MeanLogFC)) + 
    geom_density() + 
    geom_vline(xintercept = mean(this_svcaller_palmtrees$logFC), color="red") + 
    ggtitle(paste0("Test for Copy Number Increase within Palm Tree Regions (", svcallers[svi], ")")) +
    ylab("Density") +
    xlab(paste0("Mean logFC \n Random Intervals (black) vs. Palm Trees (red) \n", sprintf("p=%.2f", 1-myECDF(mean(this_svcaller_palmtrees$logFC))))) +
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToCopyNumber/CopyNumber_", svcallers[svi], "_PermutationTest.pdf"), height=4, width=6)
  
  ## Plotting
  distr = rbind(data.frame(logFC = this_svcaller_RandomDraws$logFC,
                           isPerm = "Random Intervals"),
                data.frame(logFC = this_svcaller_palmtrees$logFC,
                           isPerm = "Palm Trees"))
  
  ggplot(data=distr, aes(x=logFC, color=isPerm)) + 
    geom_density(bw=0.02, size=0.3) + 
    scale_color_manual(values = c("black", "red")) + 
    xlab("logFC") + 
    ylab("Density") + 
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          legend.title = element_blank(),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
    ggtitle(paste0("Copy Number in Palm Trees vs. Random Intervals (", svcallers[svi], ")")) +
    ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/RelationToCopyNumber/CopyNumber_", svcallers[svi], "_PalmTreesVSRandomIntervals_NotCut.pdf"),height=4, width=6)
}



# # Circos 
 # library(circlize)
# palmtrees_onlymatched$SampleXSVCaller = paste0(palmtrees_onlymatched$Sample, palmtrees_onlymatched$SVCaller)
# samplesXsvcallers = unique(palmtrees_onlymatched$SampleXSVCaller)
# for (i in 1:length(samplesXsvcallers)){
#   mycols = rand_color(n=100, luminosity = "dark")
#   
#   this_sample_this_svcaller = palmtrees_onlymatched %>% filter(SampleXSVCaller == samplesXsvcallers[i])
#   this_sample = this_sample_this_svcaller[[1,"Sample"]]
#   this_svcaller = this_sample_this_svcaller[[1,"SVCaller"]]
# 
#   pdf(file=paste0("~/Desktop/PalmTrees/Results/Figures/CircosCNVPalmTrees/", as.character(samplesXsvcallers[i]), ".pdf"))
#   circos.clear()
#   circos.par("start.degree" = 90)
#   circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
#   text(0, 0, "", cex = 1)
# 
#   ascat_cnv$logFC = log10((ascat_cnv$TumorTotalCopyNumber + 0.001) / (ascat_cnv$NormalTotalCopyNumber + 0.001))
#   ascat_cnv_bed = ascat_cnv %>% filter(Sample == this_sample) %>% dplyr::select(Chr, Start, End, logFC)
# 
#   circos.genomicTrack(as.data.frame(ascat_cnv_bed), 
#                       panel.fun = function(region, value, ...) {
#                         circos.genomicRect(region, value[[1]], ytop.column = 1, ybottom = 0, 
#                                            col = ifelse(value[[1]] > 0, "steelblue", ifelse(value[[1]] < 0, "red", "gray50")), ...)
#                       })
#   
#   
#   
#   palmtrees_of_interest = this_sample_this_svcaller
#   palmtrees_bed = palmtrees_of_interest[,c("Chr", "FirstElement", "LastElement")]
#   
#   # this is just for plotting
#   diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
#   if (sum(diff<5000000)>0){
#     palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff)/2
#     palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff)/2
#   }
#   
#   if (nrow(palmtrees_bed)>0){
#     palmtrees_bed$value1 = 1:nrow(palmtrees_of_interest)
#     palmtrees_bed$value2 = palmtrees_of_interest$PalmTreeID
#   }else{
#     palmtrees_bed$value1 = data.frame()
#     palmtrees_bed$value2 = data.frame()
#   }
#   colnames(palmtrees_bed) = c("chr", "start", "end", "value1", "value2")
#   
#   circos.genomicTrack(as.data.frame(palmtrees_bed),
#                       panel.fun = function(region, value, ...){
#                         circos.genomicRect(region, value, col=mycols[value[[1]]], border=NA, ...)
#                       },
#                       ylim=c(0,0.2), track.height=0.33*circos.par("track.height"))
# 
#   title(samplesXsvcallers[i])
#     
#   dev.off()
# 
# }
# 
