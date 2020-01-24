rm(list=ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Repeats.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")

tx_original$NRepeatsA = NA
tx_original$NRepeatsB = NA
tx_original$RepeatA = NA
tx_original$RepeatClassA = NA
tx_original$RepeatB = NA
tx_original$RepeatClassB = NA

for (i in 1:nrow(tx_original)){

  a = repeats %>% filter(hasOverlap_withChr(tx_original[[i,"ChrA"]], tx_original[[i,"PosA"]], tx_original[[i,"PosA"]],
                                            Chr, Start, End))
  if (nrow(a)>0){
    tx_original[[i,"NRepeatsA"]] = nrow(a)
    tx_original[[i,"RepeatA"]] = as.character(a[[1,"Repeat"]])
    tx_original[[i,"RepeatClassA"]] = as.character(a[[1,"Class"]])
  }
  b = repeats %>% filter(hasOverlap_withChr(tx_original[[i,"ChrB"]], tx_original[[i,"PosB"]], tx_original[[i,"PosB"]],
                                            Chr, Start, End))
  if (nrow(b)>0){
    tx_original[[i,"NRepeatsB"]] = nrow(b)
    tx_original[[i,"RepeatB"]] = as.character(b[[1,"Repeat"]])
    tx_original[[i,"RepeatClassB"]] = as.character(b[[1,"Class"]])
  }
}

tx_repeats = tx_original
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/TxVsRepeats.Rdata")
