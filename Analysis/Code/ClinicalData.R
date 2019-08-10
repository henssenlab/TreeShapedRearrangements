rm(list=ls())
library(ggplot2)
library(dplyr)
library(survminer)
library(survival)

# Open clinical_data
clinical_data = read.table("~/Desktop/PalmTrees/Data/ClinicalData.csv", header=T, sep=",")

# Prognostic relevance of palm tree (merged rearrangement calls) for all risk groups
surv_object <- Surv(time = clinical_data$Survival, event=clinical_data$hasDied)
fit <- survfit(surv_object ~ hasUnionPalmTree, data = clinical_data)
ggsurvplot(fit, data = clinical_data, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", legend.labs = c("NoPalmTrees", "PalmTrees"), ylab="Overall Survival", xlab = "Days", title = paste0("All Risk Group"))

# Prognostic relevance of palm tree (merged rearrangement calls) 
# that overlaps MYCN +- 1Mb (chr2:15080683-17087129) for MYCN-amplified NB only
mycn_subgroup = clinical_data %>% filter(Risk == "MNA")
surv_object <- Surv(time = mycn_subgroup$Survival, event=mycn_subgroup$hasDied)
fit <- survfit(surv_object ~ hasUnionPalmTreeWithMYCN, data = mycn_subgroup)
ggsurvplot(fit, data = mycn_subgroup, pval = TRUE, palette="Set1", risk.table = T, legend.title = "", legend.labs = c("NoPalmTreeAroundMYCN", "PalmTreeAroundMYCN"), ylab="Overall Survival", xlab = "Days", title = paste0("MYCN-amplified NB"))
