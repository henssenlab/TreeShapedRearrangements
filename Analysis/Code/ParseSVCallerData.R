library(dplyr)
library(tidyr)

parse_smufin = function (smufinoutput_fname){
  smufinoutput = read.table(smufinoutput_fname, header=F, sep="\t")
  d = smufinoutput[,0]
  d$Sample = lapply(as.character(smufinoutput$V9), function (s) s %>% gsub("^.*sample_","",.) %>% gsub("_intraREC.*", "", .))
  d$Cohort = NA
  d[nchar(d$Sample) == 6, "Cohort"] = "Berlin" 
  d[nchar(d$Sample) == 5, "Cohort"] = "Peifer"
  d$Sample = as.character(d$Sample)
  d$SVCaller = "Smufin"
  d$ChrA = paste0("chr", as.character(smufinoutput$V3))
  d$PosA = as.numeric(smufinoutput$V4)
  d$ChrB = paste0("chr", as.character(smufinoutput$V5))
  d$PosB = as.numeric(smufinoutput$V6)
  d
}

parse_delly = function (dellyoutput_fname, cohort_name){
  dellyoutput = read.table(dellyoutput_fname, header=F, sep="\t")
  dellyoutput = dellyoutput %>% filter(V6 == "<TRA>")
  d = dellyoutput[,0]
  d$Cohort = cohort_name
  d$SVCaller = "Delly"
  d$Sample = as.character(lapply(as.character(dellyoutput$V1), function (s) gsub("_.*", "", s)))
  d$ChrA = as.character(dellyoutput$V2)
  d$PosA = as.numeric(dellyoutput$V3)
  d$ChrB = as.character(lapply(as.character(dellyoutput$V9), function (s) s %>% gsub("^.*;CHR2=","",.) %>% gsub(";.*", "", .) %>% paste0("chr", .)))
  d$PosB = as.numeric(lapply(as.character(dellyoutput$V9), function (s) s %>% gsub("^.*;END=","",.) %>% gsub(";.*", "", .)))
  d
}

parse_novobreak = function (novobreakoutput_fname, cohort_name){
  novobreakoutput = read.table(novobreakoutput_fname, header=T, sep="\t")
  novobreakoutput = novobreakoutput %>% filter(ALT=="<TRA>")
  d = novobreakoutput[,0]
  d$Cohort = cohort_name
  d$SVCaller = "Novobreak"
  d$Sample = as.character(lapply(as.character(novobreakoutput$SAMPLE), function (s) s %>% gsub("CB", "NB", .)))
  d$ChrA = as.character(novobreakoutput$CHR.A)
  d$PosA = as.numeric(novobreakoutput$POS.A)
  d$ChrB = as.character(novobreakoutput$CHR.B)
  d$PosB = as.numeric(novobreakoutput$POS.B)
  d
}

parse_merged = function(fname){
  merged = read.table(fname, header = F, sep="\t")
  d = merged[,0]
  d$Sample = as.character(lapply(as.character(merged$V1), function (s) gsub("CB", "NB", s)))
  d$Cohort = ifelse(grepl("NB2", d$Sample), "Berlin", "Peifer")
  d$ChrA = paste0("chr", as.character(merged$V2))
  d$PosA = as.numeric(merged$V3)
  d$ChrB = paste0("chr", as.character(merged$V4))
  d$PosB = as.numeric(merged$V5)
  d$Delly = grepl("DELLY2=PASS", merged$V6)
  d$DellyLowQual = grepl("DELLY=LowQual", merged$V6)
  d$Smufin = grepl("SMUFIN=PASS", merged$V6)
  d$Svaba = grepl("SVABA=PASS", merged$V6)
  d$Novobreak = grepl("NOVOBREAK=PASS", merged$V6)
  d$Brass = grepl("BRASS=PASS", merged$V6)
  d
}

parse_pedpancan = function(fname){
  d = read.table(fname, header = T, sep=";")
  d = d %>% filter(Type == "TRA") %>% dplyr::select(-DiscordantMates, -JunctionSeqRes)
  d$SVCaller = "PedPanCanPipeline"
  d$Cohort = "PedPanCan"
  d = d %>% filter(isdefchrom2(ChrA), isdefchrom(ChrB)) 
  d$ChrA = as.character(d$ChrA)
  d$ChrB = as.character(d$ChrB)
  d
}

parse_pedpancan_allsv = function(fname){
  d = read.table(fname, header = T, sep=";")
  #d = d %>% dplyr::select(DiscordantMates, JunctionSeqRes)
  d$SVCaller = "PedPanCanPipeline"
  d$Cohort = "PedPanCan"
  d = d %>% filter(isdefchrom2(ChrA), isdefchrom(ChrB)) 
  d$ChrA = as.character(d$ChrA)
  d$ChrB = as.character(d$ChrB)
  d
}


# defines which chromosomes are taken into account, may be used to exclude e.g. chrM
isdefchrom = function(v){
  v == "chr1" | v == "chr2" | v == "chr3" | v == "chr4" | v == "chr5" |  v == "chr6" |  v == "chr7" | v == "chr8" | v == "chr9" | v == "chr10" | v == "chr11" | v == "chr12" | v == "chr13" | v == "chr14" | v == "chr15" | v == "chr16" | v == "chr17" | v == "chr18" | v == "chr19" | v == "chr20" | v == "chr21" | v == "chr22" | v == "chrX" |v == "chrY" 
}

isdefchrom2 = function(v){
  v %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
}

