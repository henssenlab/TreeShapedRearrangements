rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)
library(circlize)
setwd("~/Desktop/PalmTrees")


# Load palm tree coordinates and rearrangements
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")


# Function to create a circos plot for one sample
plot_onesample_onesvcaller = function(palmtrees, txdouble_ptinfo, txdouble, sample, svcaller, fname){

  
  # Create a random color palette
  mycols = rand_color(n=100, luminosity = "dark")

  
  # Palm trees for this sample
  palmtrees_of_interest = palmtrees %>% filter(Sample == sample, SVCaller == svcaller) %>% as.data.frame()
  
  
  # Initialize Circos plot
  pdf(file=fname)
  circos.clear()
  circos.par("start.degree" = 90)
  #circos.initializeWithIdeogram(species = "hg19")
  circos.initializeWithIdeogram(species = "hg19", plotType = c("axis", "labels"))
  text(0, 0, "", cex = 1) 
  
  
  # Create bed-type representation of palm tree regions
  palmtrees_bed = palmtrees_of_interest[,c("Chr", "FirstElement", "LastElement")] %>% as.data.frame()
  
  
  # If a palm tree region is too small to be printed in the plot, 
  # make it a little larger such that is becomes visible as a line
  diff = palmtrees_bed$LastElement-palmtrees_bed$FirstElement
  if (sum(diff<5000000)>0){
    palmtrees_bed[diff<5000000,"FirstElement"] = palmtrees_bed[diff<5000000,"FirstElement"] - (5000000-diff[diff<5000000])/2
    palmtrees_bed[diff<5000000,"LastElement"] = palmtrees_bed[diff<5000000,"LastElement"] + (5000000-diff[diff<5000000])/2
  }
  
  
  # Reformat palm tree region information to make it work with circlize grammar
  if (nrow(palmtrees_bed)>0){
    palmtrees_bed$value1 = 1:nrow(palmtrees_of_interest)
    palmtrees_bed$value2 = palmtrees_of_interest$PalmTreeID
  }else{
    palmtrees_bed$value1 = data.frame()
    palmtrees_bed$value2 = data.frame()
  }
  colnames(palmtrees_bed) = c("chr", "start", "end", "value1", "value2")
  bed_palmtree_A = data.frame()
  bed_palmtree_B = data.frame()
  if (nrow(palmtrees_bed) > 0){
    for (i in 1:nrow(palmtrees_bed)){
      palmtree_tx = txdouble_ptinfo %>% filter(PalmTreeID == as.character(palmtrees_bed[i, "value2"]))
      palmtree_tx$PalmTreeIndex = i
      palmtree_tx = palmtree_tx %>% dplyr::select(PalmTreeChrom, PalmTreePos, TargetChrom, TargetPos, PalmTreeIndex) %>% distinct()
      temp = palmtree_tx[,c("PalmTreeChrom", "PalmTreePos", "PalmTreePos", "PalmTreeIndex")]
      colnames(temp) = c("chr", "start", "end", "value1")
      bed_palmtree_A = rbind(bed_palmtree_A, temp)
      colnames(bed_palmtree_A) = c("chr", "start", "end", "value1")
      temp = palmtree_tx[,c("TargetChrom", "TargetPos", "TargetPos", "PalmTreeIndex")]
      colnames(temp) = c("chr", "start", "end", "value1")
      bed_palmtree_B = rbind(bed_palmtree_B, temp)
      colnames(bed_palmtree_B) = c("chr", "start", "end", "value1")
    }
  }
  bed_palmtree_A = as.data.frame(bed_palmtree_A)
  bed_palmtree_B = as.data.frame(bed_palmtree_B)
  
  
  # Plot palm tree regions, assign random colors to each palm tree
  circos.genomicTrack(as.data.frame(palmtrees_bed),
                      panel.fun = function(region, value, ...){
                        circos.genomicRect(region, value, col=mycols[value[[1]]], border=NA, ...)
                      },
                      ylim=c(0,0.2), track.height=0.33*circos.par("track.height"))
  
  
  # Plot all non palm tree-associated rearrangements
  nopalmtree_tx = tx_original %>% filter(Sample == sample, SVCaller == svcaller, isdefchrom(ChrA), isdefchrom(ChrB))
  bed_nopalmtree_A = nopalmtree_tx[,c("ChrA", "PosA", "PosA")]  %>% as.data.frame()
  if (nrow(bed_nopalmtree_A) > 0){
    bed_nopalmtree_A$value1 = 0
    colnames(bed_nopalmtree_A) = c("chr", "start", "end", "value1")
  }
  bed_nopalmtree_B = nopalmtree_tx[,c("ChrB", "PosB", "PosB")] %>% as.data.frame()
  if (nrow(bed_nopalmtree_B) > 0){
    bed_nopalmtree_B$value1 = 0
    colnames(bed_nopalmtree_B) = c("chr", "start", "end", "value1")
  }
  circos.genomicLink(bed_nopalmtree_A, bed_nopalmtree_B, col="lightgray")
  
  
  # Plot all palm tree-associated rearrangements in the color of the respective palm tree
  if (nrow(palmtrees_bed) > 0) circos.genomicLink(bed_palmtree_A, bed_palmtree_B, col=mycols[bed_palmtree_A[,"value1"]])
  
  
  # Make title and finish plot
  title(paste(sample, svcaller))
  dev.off()
  circos.clear()
}


# Plot circos plot for all samples and the different SV call sets
samples = unique(as.character(tx$Sample))
svcallers = unique(as.character(tx$SVCaller))
for (sai in 1:length(samples)){
  print(samples[sai])
  for (svi in 1:length(svcallers)){
    if (nrow(tx %>% filter(Sample == samples[sai], SVCaller == svcallers[svi]))>0){
      plot_onesample_onesvcaller(palmtrees, txptinfo, txdouble, samples[sai], svcallers[svi], paste0("~/Desktop/PalmTrees/Results/Figures/CircosPlots/", as.character(svcallers[svi]), "/", as.character(samples[sai]), as.character(svcallers[svi]), "PalmTreeCircos.pdf"))
    }
  }
}
