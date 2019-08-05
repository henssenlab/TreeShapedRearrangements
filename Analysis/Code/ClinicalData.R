rm(list=ls())
library(ggplot2)
library(dplyr)
library(survminer)
library(survival)
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
#chr2:16080683-16087129 MYCN hg19

svcallers = unique(tx_original$SVCaller)
svcallers = "Union"
allcallers = list()
for (i in 1:length(svcallers)){
  svcaller = svcallers[i]
  data = data.frame(Sample = unique(tx_original$Sample))
  data$Translocations = lapply(data$Sample, function (sample) tx_original %>% filter(SVCaller == svcallers[i], Sample == sample) %>% summarise(ntx=n_distinct(BPID)) %>% .$ntx) %>% as.numeric()
  #data$InterchromosomalTranslocations = lapply(data$Sample, function (sample) tx_original %>% filter(SVCaller == svcallers[i], Sample == sample, ChrA != ChrB) %>% n_distinct()) %>% as.numeric()
  #data$NumberOfPalmTrees = lapply(data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample) %>% n_distinct()) %>% as.numeric()
  data$PalmTreeTranslocations = lapply(data$Sample, function (sample) txptinfo %>% filter(SVCaller == svcallers[i], Sample == sample) %>% summarise(ntx=n_distinct(BPID)) %>% .$ntx) %>% as.numeric()
  #data$PalmTreeTranslocationsRatio = data$PalmTreeTranslocations / data$Translocations
  #data$PalmTreeInterchromosomalTranslocations = lapply(data$Sample, function (sample) ((txptinfo %>% filter(SVCaller == svcallers[i], Sample == sample, PalmTreeChrom != TargetChrom) %>% n_distinct()))) %>% as.numeric()
  #data$PalmTreeInterchrTranslocationsRatio = data$PalmTreeInterchrTranslocations / data$InterchrTranslocations
  
  write.table(data, file = paste0("~/Desktop/PalmTrees/Results/OverviewData/", svcallers[i], "_TranslocationsInPalmTrees.tsv"), quote=F, sep="\t", row.names=F)
  
  # f=ggplot(data=data %>% filter(PalmTreeTranslocations > 0), aes(x=PalmTreeTranslocationsRatio)) +
  #   geom_histogram(color=NA, fill="steelblue") + 
  #   theme_classic() +
  #   xlab("") + 
  #   scale_x_continuous(limits = c(0,1)) + 
  #   ggtitle(paste0("Proportion of Translocations within Palm Trees (", svcallers[i], ")"))
  # ggsave(filename=paste0("~/Desktop/PalmTrees/Results/OverviewData/", svcallers[i], "/TranslocationsInPalmTreesRatio.pdf"), plot=f)
  # 
  # f=ggplot(data=data %>% filter(PalmTreeTranslocations > 0), aes(x=PalmTreeInterchrTranslocationsRatio)) +
  #   geom_histogram(color=NA, fill="steelblue") + 
  #   theme_classic() +
  #   xlab("") + 
  #   scale_x_continuous(limits = c(0,1)) + 
  #   ggtitle(paste0("Proportion of Interchromosomal Translocations within Palm Trees (", svcallers[i], ")"))
  # ggsave(filename=paste0("~/Desktop/PalmTrees/Results/OverviewData/", svcallers[i], "/InterchrTranslocationsInPalmTreesRatio.pdf"), plot=f)
  # 
  data$SVCaller = svcallers[i]
  allcallers[[i]] = data
}
allcallers = do.call(rbind, allcallers)

ggplot(data=allcallers, aes(x = SVCaller, y = Translocations, fill=SVCaller)) +
  geom_violin() +
  scale_fill_brewer(palette="Set1") +
  theme_classic() + 
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 10), 
        axis.title =  element_text(face="bold", size = 10), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  scale_y_log10()

## Clinical Data

clinical_berlin = read.table("~/Desktop/PalmTrees/Data/Berlin_cohort_annotation_for_partners.csv", header=T, sep=",")
clinical_berlin$Sample = as.character(lapply(clinical_berlin$PAT_ID_BLN, function (s) gsub("CB", "NB", s)))
clinical_berlin$MNA = (clinical_berlin$GROUP %in% c("MNA", "TERT_MNA"))
clinical_berlin$Risk = clinical_berlin$GROUP
clinical_berlin$Cohort = "Berlin"
clinical_berlin$Survival = clinical_berlin$UEBDAU
clinical_berlin$EventFreeSurvival = clinical_berlin$REMDAU
clinical_berlin$hasDied = as.numeric(clinical_berlin$STATUS) == 1 # check again
clinical_berlin$hadEvent = as.numeric(clinical_berlin$STATUS) != 0 # check again
clinical_berlin = clinical_berlin %>% dplyr::select(Cohort, Sample, Risk, MNA, Survival, EventFreeSurvival, hasDied, hadEvent)
clinical_peifer = read.table("~/Desktop/PalmTrees/Data/Peifer_ClinicalData_2016_09_09_FH.csv", header=T, sep=",")
clinical_peifer$Sample = clinical_peifer$NBL.ID
clinical_peifer$MNA = (clinical_peifer$MYCN.cat %in% c("MNA_HET", "MNA"))
clinical_peifer$Risk = clinical_peifer$SUBGROUP
clinical_peifer$Cohort = "Peifer"
clinical_peifer$Survival = as.numeric(clinical_peifer$UEBDAU)
clinical_peifer$EventFreeSurvival = as.numeric(clinical_peifer$REMDAU)
clinical_peifer$hasDied = as.numeric(clinical_peifer$STATUS) == 1 # check again
clinical_peifer$hadEvent = as.numeric(clinical_peifer$STATUS) != 0 # check again
clinical_peifer = clinical_peifer %>% dplyr::select(Cohort, Sample, Risk, MNA, Survival, EventFreeSurvival, hasDied, hadEvent)

###
clinical_data = rbind(clinical_berlin, clinical_peifer)
clinical_data$isHR = (clinical_data$Risk == "MNA") | (clinical_data$Risk == "HR") | (clinical_data$Risk == "TERT_MNA") | (clinical_data$Risk == "HR_nMNA") | (clinical_data$Risk == "TERT")
clinical_data$BrassPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Brass", Sample == sample) %>% n_distinct()) %>% as.numeric()
clinical_data$SmufinPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Smufin", Sample == sample) %>% n_distinct()) %>% as.numeric()
clinical_data$NovobreakPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Novobreak", Sample == sample) %>% n_distinct()) %>% as.numeric()
clinical_data$DellyPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Delly", Sample == sample) %>% n_distinct()) %>% as.numeric()
clinical_data$SvabaPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Svaba", Sample == sample) %>% n_distinct()) %>% as.numeric()
clinical_data$AtLeastTwoCallersPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "AtLeastTwo", Sample == sample) %>% n_distinct()) %>% as.numeric()
clinical_data$UnionPalmTrees = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Union", Sample == sample) %>% n_distinct()) %>% as.numeric()

clinical_data$hasUnionPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Union", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15080683,17087129)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasUnionPalmTreeWithMYCN = clinical_data$hasUnionPalmTreeWithMYCN > 0

clinical_data$hasAtLeastTwoPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "AtLeastTwo", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasAtLeastTwoPalmTreeWithMYCN = clinical_data$hasAtLeastTwoPalmTreeWithMYCN > 0

clinical_data$hasSmufinPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Smufin", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasSmufinPalmTreeWithMYCN = clinical_data$hasSmufinPalmTreeWithMYCN > 0

clinical_data$hasDellyPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Delly", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasDellyPalmTreeWithMYCN = clinical_data$hasDellyPalmTreeWithMYCN > 0

clinical_data$hasNovobreakPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Novobreak", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasNovobreakPalmTreeWithMYCN = clinical_data$hasNovobreakPalmTreeWithMYCN > 0

clinical_data$hasBrassPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Brass", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasBrassPalmTreeWithMYCN = clinical_data$hasBrassPalmTreeWithMYCN > 0

clinical_data$hasSvabaPalmTreeWithMYCN = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Svaba", Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()
clinical_data$hasSvabaPalmTreeWithMYCN = clinical_data$hasSvabaPalmTreeWithMYCN > 0

clinical_data$hasUnionPalmTreeOnChr12 = lapply(clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == "Union", Sample == sample, Chr == "chr12") %>% n_distinct()) %>% as.numeric()
clinical_data$hasUnionPalmTreeOnChr12 = clinical_data$hasUnionPalmTreeOnChr12 > 0

clinical_data$hasPalmTree = clinical_data$SmufinPalmTrees | clinical_data$NovobreakPalmTrees | clinical_data$DellyPalmTrees | clinical_data$SvabaPalmTrees | clinical_data$AtLeastTwoCallersPalmTrees | clinical_data$UnionPalmTrees

clinical_data$hasUnionPalmTree = clinical_data$UnionPalmTrees > 0
clinical_data$hasAtLeastTwoPalmTree = clinical_data$AtLeastTwoCallersPalmTrees > 0
clinical_data$hasSmufinPalmTree = clinical_data$SmufinPalmTrees > 0
clinical_data$hasSvabaPalmTree = clinical_data$SvabaPalmTrees > 0
clinical_data$hasDellyPalmTree = clinical_data$DellyPalmTrees > 0
clinical_data$hasNovobreakPalmTree = clinical_data$NovobreakPalmTrees > 0
clinical_data$hasBrassPalmTree = clinical_data$BrassPalmTrees > 0

clinical_data[clinical_data$Risk == "HR","Risk"] = "HR_nMNA"
clinical_data[clinical_data$Risk == "TERT","Risk"] = "HR_nMNA"
clinical_data[clinical_data$Risk == "TERT_MNA","Risk"] = "MNA"

## Calling Data
callinginfo = read.csv("~/Desktop/PalmTrees/Data/EliasMail190125/NB_samples_RUN_summary.csv", header=T, sep=";")
callinginfo$Sample = gsub("CB", "NB", callinginfo$Sample)
callinginfo$AtLeastTwo = TRUE 
callinginfo$Union = TRUE 
callinginfo = callinginfo %>% filter(!(Sample %in% c("NBL47", "NBL53", "NBL54", "NBL61"))) # excludes sampels marked in the .xlsx file

## Filter Clinical Data
clinical_data = clinical_data %>% filter(Sample %in% callinginfo$Sample, !(Sample %in% c("NBL47", "NBL49", "NBL50")))


#############################################

#svcallers = unique(palmtrees$SVCaller)
#svcallers = c("Smufin", "Delly", "Svaba", "Brass", "AtLeastTwo", "Union")
svcallers = "Union"

byriskgroup = list()
for (i in 1:length(svcallers)){
  
  this_palmtrees = palmtrees %>% filter(SVCaller == svcallers[i])
  
  if (nrow(this_palmtrees) <= 2) next
  
  this_clinical_data = clinical_data %>% filter(!is.na(Risk), Risk != "NA", Sample %in% (callinginfo[callinginfo[,svcallers[i]], "Sample"])) %>% droplevels()
  this_clinical_data$ThisCallerHasPalmTree = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample) %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasPalmTreeLogical = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample) %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasPalmTree = ifelse(this_clinical_data$ThisCallerHasPalmTree,
                                                  "Palm Tree Positive",
                                                  "Palm Tree Negative")
  this_clinical_data$ThisCallerHasPalmTree = factor(this_clinical_data$ThisCallerHasPalmTree,
                                                    levels=c("Palm Tree Negative", "Palm Tree Positive"),
                                                    labels=c("No Palm Trees", "Palm Trees"))
  
  
  this_clinical_data$ThisCallerHasLargePalmTree = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample, n>=5) %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasLargePalmTreeLogical = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample, n>=5) %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasLargePalmTree = ifelse(this_clinical_data$ThisCallerHasLargePalmTree,
                                                    "Large Palm Tree Positive",
                                                    "Large Palm Tree Negative")
  this_clinical_data$ThisCallerHasLargePalmTree = factor(this_clinical_data$ThisCallerHasLargePalmTree,
                                                    levels=c("Large Palm Tree Negative", "Large Palm Tree Positive"),
                                                    labels=c("No Large Palm Trees", "Large Palm Trees"))
  
  
  #this_clinical_data$ThisCallerHasMYCNPalmTree = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 16080683, 16087129)) %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasMYCNPalmTree = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample, Chr == "chr2", hasOverlap(FirstElement, LastElement, 15000000, 17000000)) %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasMYCNPalmTree = ifelse(this_clinical_data$ThisCallerHasMYCNPalmTree,
                                                    "MYCN Palm Tree Positive",
                                                    "MYCN Palm Tree Negative")
  this_clinical_data$ThisCallerHasMYCNPalmTree = factor(this_clinical_data$ThisCallerHasMYCNPalmTree,
                                                    levels=c("MYCN Palm Tree Negative", "MYCN Palm Tree Positive"),
                                                    labels=c("No MYCN Palm Tree", "MYCN Palm Tree"))
  
  this_clinical_data$ThisCallerHasChr12PalmTree = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample, Chr == "chr12") %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasChr12PalmTree = ifelse(this_clinical_data$ThisCallerHasChr12PalmTree,
                                                        "Has chr12 Palm Tree",
                                                        "Has no chr12 Palm Tree")
  this_clinical_data$ThisCallerHasChr12PalmTree = factor(this_clinical_data$ThisCallerHasChr12PalmTree,
                                                        levels=c("Has no chr12 Palm Tree", "Has chr12 Palm Tree"),
                                                        labels=c("Has no chr12 Palm Tree", "Has chr12 Palm Tree"))
  
  this_clinical_data$ThisCallerHasChr17PalmTree = (lapply(this_clinical_data$Sample, function (sample) palmtrees %>% filter(SVCaller == svcallers[i], Sample == sample, Chr == "chr17") %>% n_distinct()) %>% as.numeric()) > 0
  this_clinical_data$ThisCallerHasChr17PalmTree = ifelse(this_clinical_data$ThisCallerHasChr17PalmTree,
                                                         "Has chr17 Palm Tree",
                                                         "Has no chr17 Palm Tree")
  this_clinical_data$ThisCallerHasChr17PalmTree = factor(this_clinical_data$ThisCallerHasChr17PalmTree,
                                                         levels=c("Has no chr17 Palm Tree", "Has chr17 Palm Tree"),
                                                         labels=c("Has no chr17 Palm Tree", "Has chr17 Palm Tree"))
  
  
  
  # How many samples have palm trees?
  this_clinical_data %>%
    mutate(Risk = factor(this_clinical_data$Risk, levels = c("ST4S", "LR", "IMR", "TERT_MNA", "TERT", "HR_nMNA", "MNA"))) %>%
    group_by(Risk) %>%
    summarise(PalmTreePositive = mean(ThisCallerHasPalmTreeLogical)) %>%
    ggplot(aes(x=Risk, y=100*PalmTreePositive)) + 
      geom_col(fill="steelblue", color="steelblue") + 
    scale_fill_brewer(palette="Set1") + 
    theme_classic() +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 14), 
          axis.title =  element_text(size = 14), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5))+
    scale_y_continuous(limits=c(0,100)) +
    guides(fill=F) + 
    xlab("") + 
    ylab("Tree-shaped Rearrangements\n[% Tumors]") + 
    #ggtitle(paste0("Palm Trees by Risk Groups (", svcallers[i], ")"))
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_RiskGroupDistribution.pdf"), width=4, height=4)

  # save that info for later
  byriskgroup[[i]] = this_clinical_data %>%
    mutate(Risk = factor(this_clinical_data$Risk, levels = c("ST4S", "LR", "IMR", "TERT_MNA", "TERT", "HR_nMNA", "MNA"))) %>%
    group_by(Risk) %>%
    summarise(SVCaller = svcallers[i], PalmTreePositive = mean(ThisCallerHasPalmTreeLogical))
  
  # all samples
  surv_object <- Surv(time = this_clinical_data$Survival, event=this_clinical_data$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data)
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_OS.pdf"), onefile=F, width=8, height=6)
  print(ggsurvplot(fit, data = this_clinical_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("All Risk Groups (", svcallers[i], ")")))
  dev.off()
  
  # all samples small figure for manuscript
  surv_object <- Surv(time = this_clinical_data$Survival, event=this_clinical_data$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data)
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_OS_SmallFigure.pdf"), onefile=F, width=4, height=4)
  print(ggsurvplot(fit, data = this_clinical_data, pval = TRUE, palette=c("firebrick2", "steelblue"), legend.title = "", legend.labs = c("No TSR", "TSR"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("All Tumors"),pval.coord = c(4000,0.1)))
  dev.off()
  
  surv_object <- Surv(time = this_clinical_data$EventFreeSurvival, event=this_clinical_data$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data)
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_EFS.pdf"), onefile=F, width=8, height=6)
  print(ggsurvplot(fit, data = this_clinical_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Event-free Survival", xlab = "Days",
                   title = paste0("All Risk Groups (", svcallers[i], ")")))
  dev.off()
  
  # only MNA
  surv_object <- Surv(time = this_clinical_data %>% filter(MNA) %>% .$Survival, event=this_clinical_data %>% filter(MNA) %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_MNAonly_OS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(MNA), pval = TRUE, palette="Set1", risk.table = T, legend.title = "",legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  surv_object <- Surv(time = this_clinical_data %>% filter(MNA) %>% .$EventFreeSurvival, event=this_clinical_data %>% filter(MNA) %>% .$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_MNAonly_EFS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(MNA), pval = TRUE, palette="Set1", risk.table=T, legend.title = "",legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Event-Free Survival", xlab = "Days",
                   title = paste0("MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  
  # only non_MNA
  surv_object <- Surv(time = this_clinical_data %>% filter(!MNA) %>% .$Survival, event=this_clinical_data %>% filter(!MNA) %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(!MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_nonMNAonly_OS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(!MNA), pval = TRUE, palette="Set1", risk.table=T, legend.title = "",legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("Non MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  surv_object <- Surv(time = this_clinical_data %>% filter(!MNA) %>% .$EventFreeSurvival, event=this_clinical_data %>% filter(!MNA) %>% .$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(!MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_nonMNAonly_EFS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(!MNA), pval = TRUE, palette="Set1", risk.table=T, legend.title = "",legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Event-free Survival", xlab = "Days",
                   title = paste0("Non MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  
  # only non_MNA --- LARGE PALMTree
  if (length(unique(this_clinical_data$ThisCallerHasLargePalmTree))>1){
    
  surv_object <- Surv(time = this_clinical_data %>% filter(!MNA) %>% .$Survival, event=this_clinical_data %>% filter(!MNA) %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasLargePalmTree, data = this_clinical_data %>% filter(!MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_LargePT_nonMNAonly_OS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(!MNA), pval = TRUE, palette="Set1", risk.table=T, legend.title = "",legend.labs = c("NoLargePalmTrees", "LargePalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("Non MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  surv_object <- Surv(time = this_clinical_data %>% filter(!MNA) %>% .$EventFreeSurvival, event=this_clinical_data %>% filter(!MNA) %>% .$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasLargePalmTree, data = this_clinical_data %>% filter(!MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_LargePT_nonMNAonly_EFS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(!MNA), pval = TRUE, palette="Set1", risk.table=T, legend.title = "",legend.labs = c("NoLargePalmTrees", "LargePalmTrees"), ylab="Event-free Survival", xlab = "Days",
                   title = paste0("Non MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  }
  
  # only HR
  surv_object <- Surv(time = this_clinical_data %>% filter(isHR) %>% .$Survival, event=this_clinical_data %>% filter(isHR) %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(isHR))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_HRonly_OS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(isHR), pval = TRUE, palette="Set1", risk.table=T, legend.title = "", legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("High-Risk Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  
  surv_object <- Surv(time = this_clinical_data %>% filter(isHR) %>% .$EventFreeSurvival, event=this_clinical_data %>% filter(isHR) %>% .$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(isHR))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_HRonly_EFS.pdf"), onefile=F, width=8, height=6)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(isHR), pval = TRUE, palette="Set1", risk.table=T, 
                   legend.title = "",legend.labs = c("NoPalmTrees", "PalmTrees"), 
                   ylab="Event-free Survival", xlab = "Days",
                   title = paste0("High-Risk Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  # only HR_nMNA
  surv_object <- Surv(time = this_clinical_data %>% filter(Risk == "HR_nMNA") %>% .$Survival, event=this_clinical_data %>% filter(Risk == "HR_nMNA") %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(Risk == "HR_nMNA"))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_HRnMNAonly_OS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(Risk == "HR_nMNA"), pval = TRUE, palette="Set1", risk.table=T, legend.title = "", legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("HR_nMNA (", svcallers[i], ")")))
  dev.off()
  
  
  surv_object <- Surv(time = this_clinical_data %>% filter(Risk == "HR_nMNA") %>% .$EventFreeSurvival, event=this_clinical_data %>% filter(Risk == "HR_nMNA") %>% .$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasPalmTree, data = this_clinical_data %>% filter(Risk == "HR_nMNA"))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_HR_nMNA_EFS.pdf"), onefile=F, width=8, height=6)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(Risk == "HR_nMNA"), pval = TRUE, palette="Set1", risk.table=T, 
                   legend.title = "",legend.labs = c("NoPalmTrees", "PalmTrees"), 
                   ylab="Event-free Survival", xlab = "Days",
                   title = paste0("HR_nMNA (", svcallers[i], ")")))
  dev.off()
  
  # only MNA MYCN PalmTrees
  if (length(unique(this_clinical_data$ThisCallerHasMYCNPalmTree))>1){
  surv_object <- Surv(time = this_clinical_data %>% filter(MNA) %>% .$Survival, event=this_clinical_data %>% filter(MNA) %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasMYCNPalmTree, data = this_clinical_data %>% filter(MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_MNAonly_MYCNPT_OS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(MNA), pval = TRUE, palette="Set1", risk.table = T, legend.title = "",legend.labs = c("NoMYCNPalmTrees", "MYCNPalmTrees"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  surv_object <- Surv(time = this_clinical_data %>% filter(MNA) %>% .$EventFreeSurvival, event=this_clinical_data %>% filter(MNA) %>% .$hadEvent)
  fit <- survfit(surv_object ~ ThisCallerHasMYCNPalmTree, data = this_clinical_data %>% filter(MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_MNAonly_MYCNPT_EFS.pdf"), width=8, height=6, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(MNA), pval = TRUE, palette="Set1", risk.table=T, legend.title = "",legend.labs = c("NoMYCNPalmTrees", "MYCNPalmTrees"), ylab="Event-Free Survival", xlab = "Days",
                   title = paste0("MYCN Amplified Neuroblastoma (", svcallers[i], ")")))
  dev.off()
  
  # Small Figure For Manuscript
  surv_object <- Surv(time = this_clinical_data %>% filter(MNA) %>% .$Survival, event=this_clinical_data %>% filter(MNA) %>% .$hasDied)
  fit <- survfit(surv_object ~ ThisCallerHasMYCNPalmTree, data = this_clinical_data %>% filter(MNA))
  pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_MNAonly_MYCNPT_OS_SmallFigure.pdf"), width=4, height=4, onefile=F)
  print(ggsurvplot(fit, data = this_clinical_data %>% filter(MNA), pval = TRUE, palette=c("firebrick2", "steelblue"), risk.table = F, legend.title = "",legend.labs = c("No TSR at MYCN", "TSR at MYCN"), ylab="Overall Survival", xlab = "Days",
                   title = paste0("MNA"),pval.coord = c(2200,0.1)))
  dev.off()
  }
  
  # Chr 12 Palm Trees
  if (length(unique(this_clinical_data$ThisCallerHasChr12PalmTree))>1){
    surv_object <- Surv(time = this_clinical_data %>% .$Survival, event=this_clinical_data %>% .$hasDied)
    fit <- survfit(surv_object ~ ThisCallerHasChr12PalmTree, data = this_clinical_data)
    pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_Chr12PT_OS.pdf"), width=8, height=6, onefile=F)
    print(ggsurvplot(fit, data = this_clinical_data , pval = TRUE, palette="Set1", risk.table = T, legend.title = "",legend.labs = c("No Chr12 PalmTree", "Chr12 Palm Tree"), ylab="Overall Survival", xlab = "Days",
                     title = paste0("Neuroblastoma (", svcallers[i], ")")))
    dev.off()
    
    surv_object <- Surv(time = this_clinical_data %>% .$EventFreeSurvival, event=this_clinical_data %>% .$hadEvent)
    fit <- survfit(surv_object ~ ThisCallerHasChr12PalmTree, data = this_clinical_data)
    pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_Chr12PT_EFS.pdf"), width=8, height=6, onefile=F)
    print(ggsurvplot(fit, data = this_clinical_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "",legend.labs = c("No Chr12 PalmTree", "Chr12 Palm Tree"), ylab="Event-free Survival", xlab = "Days",
                     title = paste0("Neuroblastoma (", svcallers[i], ")")))
    dev.off()
  }
  
  # Chr 17 Palm Trees
  if (length(unique(this_clinical_data$ThisCallerHasChr17PalmTree))>1){
    surv_object <- Surv(time = this_clinical_data %>% .$Survival, event=this_clinical_data %>% .$hasDied)
    fit <- survfit(surv_object ~ ThisCallerHasChr17PalmTree, data = this_clinical_data)
    pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_Chr17PT_OS.pdf"), width=8, height=6, onefile=F)
    print(ggsurvplot(fit, data = this_clinical_data , pval = TRUE, palette="Set1", risk.table = T, legend.title = "",legend.labs = c("No Chr17 PalmTree", "Chr17 Palm Tree"), ylab="Overall Survival", xlab = "Days",
                     title = paste0("Neuroblastoma (", svcallers[i], ")")))
    dev.off()
    
    surv_object <- Surv(time = this_clinical_data %>% .$EventFreeSurvival, event=this_clinical_data %>% .$hadEvent)
    fit <- survfit(surv_object ~ ThisCallerHasChr17PalmTree, data = this_clinical_data)
    pdf(paste0("~/Desktop/PalmTrees/Results/Figures/Clinical/", svcallers[i], "_Chr17PT_EFS.pdf"), width=8, height=6, onefile=F)
    print(ggsurvplot(fit, data = this_clinical_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "",legend.labs = c("No Chr17 PalmTree", "Chr17 Palm Tree"), ylab="Event-free Survival", xlab = "Days",
                     title = paste0("Neuroblastoma (", svcallers[i], ")")))
    dev.off()
  }
  
}

byriskgroupall = do.call(rbind, byriskgroup)

ggplot(data=do.call(rbind, byriskgroup), aes(x=Risk, y=100*PalmTreePositive, color=SVCaller)) + 
  geom_jitter(width = 0.3, height=0) + 
  stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", size = 0.2, color = "black", width=0.5) + 
  scale_color_brewer(palette="Set1") + 
  theme_classic() +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 10), 
        axis.title =  element_text(face="plain", size = 10), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5)) + 
  scale_y_continuous(limits=c(0,100)) +
  xlab("") + 
  ylab("Palm Tree Positive Samples [%]") + 
  ggtitle(paste0("Palm Trees by Risk Groups (All SV Callers)")) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Clinical/RecurrenceAllCallers.pdf") 

# fit1 <- survfit(surv_object ~ hasPalmTree + isHR, data = clinical_data)
# system("rm ~/Desktop/PalmTrees/Results/Figures/Clinical/KaplanMeier_PTXHR.pdf")
# pdf("~/Desktop/PalmTrees/Results/Figures/Clinical/KaplanMeier_PTXHR.pdf", onefile=F, width=20, height=10)
# ggsurvplot(fit1, data = clinical_data, palette="Set1", pval=F)
# dev.off()

######

matchnames = read.table("~/Desktop/PalmTrees/Data/BankID_fromRocioEmail180612.csv", sep=";", header=T)
matchnames = matchnames %>% mutate(Sample = gsub("CB", "NB", Sample))

clinical_data %>% mutate(Patient=Sample) %>% 
  left_join(matchnames, by="Sample") %>%
  dplyr::select(Patient, BankID,Risk, Survival, hasDied, EventFreeSurvival, hadEvent,
                hasAtLeastTwoPalmTree, hasNovobreakPalmTree, hasBrassPalmTree, hasSvabaPalmTree, hasSmufinPalmTree, hasDellyPalmTree, hasUnionPalmTree,
                hasAtLeastTwoPalmTreeWithMYCN, hasNovobreakPalmTreeWithMYCN, hasBrassPalmTreeWithMYCN, hasSvabaPalmTreeWithMYCN, hasSmufinPalmTreeWithMYCN, hasDellyPalmTreeWithMYCN, hasUnionPalmTreeWithMYCN) %>% 
  write.table("~/Desktop/PalmTrees/Results/Tables/PalmTreesClinicalData.csv", col.names=T, row.names=F, sep=";")

clinical_data %>% 
  left_join(matchnames, by="Sample") %>%
  mutate(Patient=Sample) %>% 
  filter(Risk == "MNA") %>%
  dplyr::select(Patient, BankID, Risk, Survival, hasDied, EventFreeSurvival, hadEvent, 
                hasAtLeastTwoPalmTree, hasNovobreakPalmTree, hasBrassPalmTree, hasSvabaPalmTree, hasSmufinPalmTree, hasDellyPalmTree, hasUnionPalmTree,
                hasAtLeastTwoPalmTreeWithMYCN, hasNovobreakPalmTreeWithMYCN, hasBrassPalmTreeWithMYCN, hasSvabaPalmTreeWithMYCN, hasSmufinPalmTreeWithMYCN, hasDellyPalmTreeWithMYCN, hasUnionPalmTreeWithMYCN) %>% 
  write.table("~/Desktop/PalmTrees/Results/Tables/PalmTreesClinicalData_onlyMNA.csv", col.names=T, row.names=F, sep=";")

