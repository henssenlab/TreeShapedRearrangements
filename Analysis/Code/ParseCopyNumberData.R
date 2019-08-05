library(dplyr)
library(tidyr)
library(GenomicRanges)

system("mkdir ~/Desktop/PalmTrees/Analysis/WorkspaceData/WGSCopyNumberVariants/")
system("cp ~/Desktop/PalmTrees/Data/WGSCopyNumberVariants/*.CNVs ~/Desktop/PalmTrees/Analysis/WorkspaceData/WGSCopyNumberVariants/")
system("for i in ~/Desktop/PalmTrees/Analysis/WorkspaceData/WGSCopyNumberVariants/*.CNVs; do awk '{print $0 \"\\t\" FILENAME}' $i > ${i}.tmp ; mv ${i}.tmp ${i} ; done")
system("awk FNR-1 ~/Desktop/PalmTrees/Analysis/WorkspaceData/WGSCopyNumberVariants/*.CNVs > ~/Desktop/PalmTrees/Analysis/WorkspaceData/AllSamples.CNVs")
system("rm -r ~/Desktop/PalmTrees/Analysis/WorkspaceData/WGSCopyNumberVariants/")

cnv = read.table("~/Desktop/PalmTrees/Analysis/WorkspaceData/AllSamples.CNVs", header=F, sep="\t",
                  col.names = c("Chr", "Start", "End", "CopyNumber", "CNStatus", "Sample")) %>% as_tibble()
cnv$Chr = paste0("chr",as.character(cnv$Chr))
cnv$Sample = as.character(lapply(cnv$Sample, function(s) gsub("CB", "NB", gsub("-.*", "", gsub(".*freec.", "", s)))))

cnv_gr = makeGRangesFromDataFrame(cnv,
                                  seqnames.field = "Chr", 
                                  start.field = "Start",
                                  end.field = "End",
                                  keep.extra.columns = TRUE)

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/CNVData.Rdata")

system("mkdir ~/Desktop/PalmTrees/Analysis/WorkspaceData/ALL_patients_neuroblastoma_results_ASCAT_CNV/")
system("cp ~/Desktop/PalmTrees/Data/ALL_patients_neuroblastoma_results_ASCAT_CNV/*.vcf.gz  ~/Desktop/PalmTrees/Analysis/WorkspaceData/ALL_patients_neuroblastoma_results_ASCAT_CNV/")
system("gunzip ~/Desktop/PalmTrees/Analysis/WorkspaceData/ALL_patients_neuroblastoma_results_ASCAT_CNV/*.vcf.gz")
system("for i in ~/Desktop/PalmTrees/Analysis/WorkspaceData/ALL_patients_neuroblastoma_results_ASCAT_CNV/*.vcf; do awk '{print $0 \"\\t\" FILENAME}' $i > ${i}.tmp ; mv ${i}.tmp ${i} ; done")
system("awk FNR-1 ~/Desktop/PalmTrees/Analysis/WorkspaceData/ALL_patients_neuroblastoma_results_ASCAT_CNV/*.vcf > ~/Desktop/PalmTrees/Analysis/WorkspaceData/AllSamples_ASCAT_CNV.vcf")
system("rm -r ~/Desktop/PalmTrees/Analysis/WorkspaceData/ALL_patients_neuroblastoma_results_ASCAT_CNV/")

ascat_cnv = read.table("~/Desktop/PalmTrees/Analysis/WorkspaceData/AllSamples_ASCAT_CNV.vcf", comment.char="#")
ascat_cnv$Chr = as.character(paste0("chr", ascat_cnv$V1))
ascat_cnv$Start = ascat_cnv$V2
ascat_cnv$End = as.integer(lapply(ascat_cnv$V8, function (s) gsub(".*END=", "", s)))
ascat_cnv$NormalTotalCopyNumber = as.numeric(lapply(ascat_cnv$V10, function (s) gsub(".*:", "", gsub(":[0123456789]*$","",s))))
ascat_cnv$NormalMinorAlleleCopyNumber = as.numeric(lapply(ascat_cnv$V10, function (s) gsub(".*:", "", s)))
ascat_cnv$TumorTotalCopyNumber = as.numeric(lapply(ascat_cnv$V11, function (s) gsub(".*:", "", gsub(":[0123456789]*$","",s))))
ascat_cnv$TumorMinorAlleleCopyNumber = as.numeric(lapply(ascat_cnv$V11, function (s) gsub(".*:", "", s)))
ascat_cnv$Sample = as.character(lapply(ascat_cnv$V12, function (s) gsub("CB", "NB", gsub(".*/", "", gsub(".copynumber.vcf", "", s)))))
ascat_cnv = ascat_cnv %>% as_tibble() %>% select(Sample,Chr,Start,End,NormalTotalCopyNumber,NormalMinorAlleleCopyNumber,TumorTotalCopyNumber,TumorMinorAlleleCopyNumber)

ascat_cnv_gr = makeGRangesFromDataFrame(ascat_cnv, keep.extra.columns = T)

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
