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
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")

for (i in 1:nrow(palmtrees)){
  
  overlapping_trees = palmtrees %>% 
    filter(PalmTreeID != palmtrees[[i, "PalmTreeID"]],
           Chr == palmtrees[[i, "Chr"]], 
           hasOverlap(FirstElement, LastElement, palmtrees[[i, "FirstElement"]], palmtrees[[i, "LastElement"]]))
  print("========================================================================")
  print(paste0("Which Palm Trees overlap with ", as.character(palmtrees[[i, "PalmTreeID"]]), " ?"))
  print(overlapping_trees %>% dplyr::select(PalmTreeID))
  print(paste0("Which Palm Trees overlap with ", as.character(palmtrees[[i, "PalmTreeID"]]), " IN THE SAME PATIENT?"))
  print(overlapping_trees %>% filter(Sample == palmtrees[[i, "Sample"]]) %>% dplyr::select(PalmTreeID))
  print(paste0("Which Palm Trees overlap with ", as.character(palmtrees[[i, "PalmTreeID"]]), " IN ANOTHER PATIENT?"))
  print(overlapping_trees %>% filter(Sample != palmtrees[[i, "Sample"]]) %>% dplyr::select(PalmTreeID))
  print(paste0("Which other Patients have palm tree overlaping with ", as.character(palmtrees[[i, "PalmTreeID"]])))
  print(overlapping_trees %>% filter(Sample != palmtrees[[i, "Sample"]]) %>% dplyr::select(Sample) %>% distinct())
  print("========================================================================")
  
  this_recurrent_trees = data.frame(
    PalmTreeID = palmtrees[[i, "PalmTreeID"]],
    AllOverlappingTrees = nrow(overlapping_trees),
    AllOverlappingTreesSamePatient = nrow(overlapping_trees %>% filter(Sample == palmtrees[[i, "Sample"]])),
    AllOverlappingTreesOtherPatientsAnyCaller = nrow(overlapping_trees %>% filter(Sample != palmtrees[[i, "Sample"]]) %>% dplyr::select(Sample) %>% distinct())
  )
  
  if (i==1){
    recurrent_trees = this_recurrent_trees
  }else{
    recurrent_trees = rbind(recurrent_trees, this_recurrent_trees)
  }
  
}

recurrent_trees %>% View()

f=ggplot(data=recurrent_trees, aes(x=AllOverlappingTreesOtherPatientsAnyCaller)) + 
  geom_histogram(fill="steelblue") + 
  theme_minimal() + 
  xlab("How many patients have an overlapping palm tree?") +
  ylab("Palm tree count") +
  scale_x_continuous(breaks=2*c(0:5))
ggsave("~/Desktop/PalmTrees/Results/Figures/PalmTreeRecurrence/PatientsWithOverlappingPalmTrees.pdf", f, height=4, width=4)

sprintf("Probability that one tree is also found in another SV caller: %0.1f%%", 100*mean(recurrent_trees$AllOverlappingTreesSamePatient > 0)) 
sprintf("Probability that one tree is found by all callers: %0.1f%%", 100*mean(recurrent_trees$AllOverlappingTreesSamePatient == 2)) 
sprintf("Probability that a tree overlaps with a tree in another patient: %0.1f%%", 100*mean(recurrent_trees$AllOverlappingTreesOtherPatientsAnyCaller > 0))
