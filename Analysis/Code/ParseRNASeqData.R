rm(list = ls())
library(dplyr)
library(tidyr)
berlin_rna = read.table("~/Desktop/PalmTrees/Data/berlin_cohort_rnaseq_fpkm.txt", header=TRUE)
colnames(berlin_rna) = lapply(colnames(berlin_rna), function (s) gsub("CB", "NB", s))
peifer_rna = read.table("~/Desktop/PalmTrees/Data/peifer_54nb_fpkms.txt", header=TRUE)
fpkmToTpm <- function(fpkm) exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
peifer_rna_tpm = as.data.frame(lapply(peifer_rna[,2:55], fpkmToTpm))
peifer_rna_tpm = cbind(peifer_rna$Gene, peifer_rna_tpm)
colnames(peifer_rna_tpm)[1] = "Gene"
berlin_rna_tpm = as.data.frame(lapply(berlin_rna[,2:56], fpkmToTpm))
berlin_rna_tpm = cbind(berlin_rna$Gene, berlin_rna_tpm)
colnames(berlin_rna_tpm)[1] = "Gene"
peifer_rna_tpm_long = peifer_rna_tpm %>% gather(key="Sample", value="TPM", 2:55)
peifer_rna_tpm_long$Cohort = "Peifer"
berlin_rna_tpm_long = berlin_rna_tpm %>% gather(key="Sample", value="TPM", 2:56)
berlin_rna_tpm_long$Cohort = "Berlin"
peifer_rna_tpm_means = peifer_rna_tpm_long %>% group_by(Gene, Cohort) %>% summarise(MeanTPM = mean(TPM), SdTPM=sd(TPM), MedianTPM = median(TPM), MadTPM = mad(TPM)) %>% ungroup()
berlin_rna_tpm_means = berlin_rna_tpm_long %>% group_by(Gene, Cohort) %>% summarise(MeanTPM = mean(TPM), SdTPM=sd(TPM), MedianTPM = median(TPM), MadTPM = mad(TPM)) %>% ungroup()
peifer_rna_tpm_long$log2TPM = log2(peifer_rna_tpm_long$TPM + 0.000001)
berlin_rna_tpm_long$log2TPM = log2(berlin_rna_tpm_long$TPM + 0.000001)
peifer_rna_tpm_means$log2MeanTPM = log2(peifer_rna_tpm_means$MeanTPM + 0.000001)
berlin_rna_tpm_means$log2MeanTPM = log2(berlin_rna_tpm_means$MeanTPM + 0.000001)
rna_tpm = rbind(berlin_rna_tpm_long, peifer_rna_tpm_long)
rna_tpm_means = rbind(berlin_rna_tpm_means, peifer_rna_tpm_means)

# modz = function(v){
#   (v-median(v))/(median(abs(v-median(v)))+0.000001)
# }
# rna_tpm = rna_tpm %>% group_by(Gene) %>% mutate(ModZ = modz(TPM)) %>% ungroup()

#plot(x = rna_tpm %>% filter(Gene == "MYCN") %>% .$TPM, y = rna_tpm %>% filter(Gene == "MYCN") %>% .$ModZ)
#plot(x = rna_tpm %>% filter(Gene == "GAPDH") %>% .$TPM, y = rna_tpm %>% filter(Gene == "GAPDH") %>% .$ModZ)
#plot(x = rna_tpm %>% filter(Gene == "ACTB") %>% .$TPM, y = rna_tpm %>% filter(Gene == "ACTB") %>% .$ModZ)

#rna_tpm$MeanTPM = mapply(function (this_gene, this_cohort) rna_tpm_means[rna_tpm_means$Gene == this_gene & rna_tpm_means$Cohort == this_cohort, "MeanTPM"], rna_tpm$Gene, rna_tpm$Cohort)
#rna$logFC = log10( (rna_tpm$TPM + 0.000001) / (rna_tpm$MeanTPM + 0.000001) )
rm(berlin_rna, peifer_rna, berlin_rna_tpm, peifer_rna_tpm)
#rm(list=setdiff(ls(), c("rna_tpm")))
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/RNASeqData.Rdata")
