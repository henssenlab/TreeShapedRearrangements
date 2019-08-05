# cp ~/Desktop/eccDNADataAnalysis/*Nov2018/*FilterFullCrc2x.bed ~/Desktop/PalmTrees/Data/WGSCircles/
# cd ~/Desktop/PalmTrees/Data/WGSCircles/
# for i in *FilterFullCrc2x.bed; do awk '{print $0 "\t" FILENAME}' $i > ${i}.tmp ; mv ${i}.tmp ${i} ; done
# awk FNR-1 *FilterFullCrc2x.bed > AllSamples_FilterFullCrc2x.bed
library(dplyr)
library(tidyr)
wgscircles = read.table("~/Desktop/PalmTrees/Data/WGSCircles/AllSamples_FilterFullCrc2x.bed", sep="\t", header=F)
colnames(wgscircles) = c("Chr", "Start", "End", "CountRatioTumorNormal", "Length", "CountTumor", "CountNormal", "Sample")
wgscircles$Chr = paste0("chr",as.character(wgscircles$Chr))
wgscircles$Sample = lapply(wgscircles$Sample, function (s) gsub("_.*", "", s))

