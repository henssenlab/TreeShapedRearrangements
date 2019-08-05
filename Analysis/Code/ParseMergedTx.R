rm(list=ls())
library(dplyr)
library(tidyr)

parse_caller_info = function(tx_row, SVCaller){
  if (nrow(tx_row)!=1) stop("Exactly one row must be handed over.")
  s = tx_row$SVCallerInfo
  s = strsplit(as.character(s), split=";")[[1]]
  s = s[grepl(SVCaller, s)]
  s = gsub(paste0(SVCaller, "="), "", s)
  s = data.frame(s)
  s = s %>% tidyr::separate("s", into=c("Filter", "ID", "Orientation"), sep=",")
  return(cbind(do.call(rbind, rep(list(tx_row), nrow(s))), s))
}

parse_caller_info_precision = function(tx_row, SVCaller){
  if (nrow(tx_row)!=1) stop("Exactly one row must be handed over.")
  s = tx_row$SVCallerInfo
  s = strsplit(as.character(s), split=";")[[1]]
  s = s[grepl(SVCaller, s)]
  s = gsub(paste0(SVCaller, "="), "", s)
  s = data.frame(s)
  s = s %>% tidyr::separate("s", into=c("Filter", "ID", "Orientation", "Precision"), sep=",")
  return(cbind(do.call(rbind, rep(list(tx_row), nrow(s))), s))
}


get_svcaller_tx = function(tx, SVCaller){
  this_svcaller_tx = tx %>% filter(grepl(SVCaller, SVCallerInfo))
  parsed = list()
  for (i in 1:nrow(this_svcaller_tx)){
    parsed[[i]] = parse_caller_info(this_svcaller_tx[i,], SVCaller)
  }
  parsed = do.call(rbind, parsed) %>% as_tibble()
  return(parsed)
}

get_svcaller_tx_precision = function(tx, SVCaller){
  this_svcaller_tx = tx %>% filter(grepl(SVCaller, SVCallerInfo))
  parsed = list()
  for (i in 1:nrow(this_svcaller_tx)){
    parsed[[i]] = parse_caller_info_precision(this_svcaller_tx[i,], SVCaller)
  }
  parsed = do.call(rbind, parsed) %>% as_tibble()
  return(parsed)
}

makeBPID = function(Sample, ChrA, PosA, ChrB, PosB){
  chrord = c("1","2","3","4","5","6","7","8","9","10","11",
             "12","13","14","15","16","17","18","19","20",
             "21","22","X","Y","M") %>% paste0("chr",.)
  if (!(ChrA %in% chrord)) stop("ChrA not a chromsome")
  if (!(ChrB %in% chrord)) stop("ChrB not a chromsome")
  if (is.na(ChrA) | is.na(PosA) | is.na(ChrB) | is.na(PosB)) stop("NA not allowed as function input")
  if (ChrA == ChrB){
    return(ifelse(PosA <= PosB, paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB), paste0(Sample, "_", ChrB, ":", PosB, "_", ChrA, ":", PosA)))
  } else {
    return(ifelse(which(ChrA == chrord) < which(ChrB == chrord), paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB), paste0(Sample, "_", ChrB, ":", PosB, "_", ChrA, ":", PosA)))
  }
}

tx_circleseq = read.table("~/Desktop/PalmTrees/Data/EliasMail190707/ALL_PATIENTS_allCIRCLES_circleseq_resultscirculome_CLASSIFICATION_TRANS_circletype_strandinfo_precisioninfo_clinicaldatasubset_7.7.19",
                          sep="\t",
                          header=F) %>% as_tibble()
colnames(tx_circleseq) = c("Sample", "ChrA", "PosA", "ChrB", "PosB", "SVCallerInfo", "MergeInfo", "CircleGenomeClass", "CircleWindow", 
                           "NumCirclesAB", "eccDNA_NumCirclesAB", "ecDNA_NumCirclesAB")
tx_circleseq = tx_circleseq %>%
  mutate(CircleMethod = "CircleSeq",
         Sample = gsub("CB", "NB", Sample),
         Cohort = ifelse(grepl("NBL", Sample), "Peifer", "Berlin"),
         ChrA = paste0("chr", ChrA),
         PosA = as.integer(PosA),
         ChrB = paste0("chr", ChrB), 
         PosB = as.integer(PosB),
         NumCirclesAB = ifelse(NumCirclesAB == "None", "0-0", as.character(NumCirclesAB)),
         eccDNA_NumCirclesAB = ifelse(eccDNA_NumCirclesAB == "None", "0-0", gsub("eccDNA:", "", as.character(eccDNA_NumCirclesAB))),
         ecDNA_NumCirclesAB = ifelse(ecDNA_NumCirclesAB == "None", "0-0", gsub("ecDNA:", "", as.character(ecDNA_NumCirclesAB)))) %>%
  separate(NumCirclesAB, into=c("NumCirclesA", "NumCirclesB"), sep="-") %>%
  separate(eccDNA_NumCirclesAB, into=c("eccDNA_NumCirclesA", "eccDNA_NumCirclesB"), sep="-") %>% 
  separate(ecDNA_NumCirclesAB, into=c("ecDNA_NumCirclesA", "ecDNA_NumCirclesB"), sep="-") %>%
  #mutate(CircleWindow = ifelse(CircleWindow == "None", NA, as.character(CircleWindow))) %>%
  #separate(CircleWindow, into=c("CircleWindowChr", "CircleWindowStart", "CircleWindowEnd"), sep="[-:]") %>%
  separate(MergeInfo, into=c("MergeInfo_Filter", "MergeInfo_nbkps", "MergeInfo_ReadSupport"), sep=",") %>%
  mutate(Svaba = grepl("SVABA", SVCallerInfo),
         Smufin = grepl("SMUFIN", SVCallerInfo),
         Delly = grepl("DELLY2", SVCallerInfo),
         Brass = grepl("BRASS", SVCallerInfo),
         Novobreak = grepl("NOVOBREAK", SVCallerInfo)) %>%
  filter(Sample != "NBL50", Sample != "NBL49", Sample != "NBL47") %>%
  filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
tx_circleseq$BPID = mapply(makeBPID, tx_circleseq$Sample, tx_circleseq$ChrA, tx_circleseq$PosA, tx_circleseq$ChrB, tx_circleseq$PosB)

delly_circleseq = get_svcaller_tx_precision(tx_circleseq, "DELLY2") %>% mutate(SVCaller = "Delly")
svaba_circleseq = get_svcaller_tx_precision(tx_circleseq, "SVABA") %>% mutate(SVCaller = "Svaba")
smufin_circleseq = get_svcaller_tx_precision(tx_circleseq, "SMUFIN") %>% mutate(SVCaller = "Smufin")
novobreak_circleseq = get_svcaller_tx_precision(tx_circleseq, "NOVOBREAK") %>% mutate(SVCaller = "Novobreak")
brass_circleseq = get_svcaller_tx_precision(tx_circleseq, "BRASS") %>% mutate(SVCaller = "Brass")
union_circleseq = tx_circleseq %>% mutate(SVCaller = "Union")
atleasttwo_circleseq = tx_circleseq %>% filter((Smufin+Delly+Svaba+Brass+Novobreak)>=2)  %>% mutate(SVCaller = "AtLeastTwo")
tx_circleseq = bind_rows(
  delly_circleseq %>% dplyr::select(-Filter, -ID, -Orientation),
  smufin_circleseq %>% dplyr::select(-Filter, -ID, -Orientation),
  svaba_circleseq %>% dplyr::select(-Filter, -ID, -Orientation),
  novobreak_circleseq %>% dplyr::select(-Filter, -ID, -Orientation),
  brass_circleseq %>% dplyr::select(-Filter, -ID, -Orientation),
  union_circleseq,
  atleasttwo_circleseq
)

### WGS ###

tx_wgscircles = read.table("~/Desktop/PalmTrees/Data/EliasMail190707/ALL_PATIENTS_allCIRCLES_WGS_resultscirculome_CLASSIFICATION_TRANS_circletype_strandinfo_precisioninfo_clinicaldatasubset_7.7.19",
                           sep="\t",
                           header=F)
colnames(tx_wgscircles) = c("Sample", "ChrA", "PosA", "ChrB", "PosB", "SVCallerInfo", "MergeInfo", "CircleGenomeClass", "CircleWindow", 
                           "NumCirclesAB", "eccDNA_NumCirclesAB", "ecDNA_NumCirclesAB")
tx_wgscircles = tx_wgscircles %>%
  mutate(CircleMethod = "WGSCircles",
         Sample = gsub("CB", "NB", Sample),
         Cohort = ifelse(grepl("NBL", Sample), "Peifer", "Berlin"),
         ChrA = paste0("chr", ChrA),
         PosA = as.integer(PosA),
         ChrB = paste0("chr", ChrB), 
         PosB = as.integer(PosB),
         NumCirclesAB = ifelse(NumCirclesAB == "None", "0-0", as.character(NumCirclesAB)),
         eccDNA_NumCirclesAB = ifelse(eccDNA_NumCirclesAB == "None", "0-0", gsub("eccDNA:", "", as.character(eccDNA_NumCirclesAB))),
         ecDNA_NumCirclesAB = ifelse(ecDNA_NumCirclesAB == "None", "0-0", gsub("ecDNA:", "", as.character(ecDNA_NumCirclesAB)))) %>%
  separate(NumCirclesAB, into=c("NumCirclesA", "NumCirclesB"), sep="-") %>%
  separate(eccDNA_NumCirclesAB, into=c("eccDNA_NumCirclesA", "eccDNA_NumCirclesB"), sep="-") %>% 
  separate(ecDNA_NumCirclesAB, into=c("ecDNA_NumCirclesA", "ecDNA_NumCirclesB"), sep="-") %>%
  #mutate(CircleWindow = ifelse(CircleWindow == "None", NA, as.character(CircleWindow))) %>%
  #separate(CircleWindow, into=c("CircleWindowChr", "CircleWindowStart", "CircleWindowEnd"), sep="[-:]") %>%
  separate(MergeInfo, into=c("MergeInfo_Filter", "MergeInfo_nbkps", "MergeInfo_ReadSupport"), sep=",") %>%
  mutate(Svaba = grepl("SVABA", SVCallerInfo),
         Smufin = grepl("SMUFIN", SVCallerInfo),
         Delly = grepl("DELLY2", SVCallerInfo),
         Brass = grepl("BRASS", SVCallerInfo),
         Novobreak = grepl("NOVOBREAK", SVCallerInfo)) %>% 
  filter(Sample != "NBL50", Sample != "NBL49", Sample != "NBL47") %>%
  filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
tx_wgscircles$BPID = mapply(makeBPID, tx_wgscircles$Sample, tx_wgscircles$ChrA, tx_wgscircles$PosA, tx_wgscircles$ChrB, tx_wgscircles$PosB)

delly_wgscircles = get_svcaller_tx_precision(tx_wgscircles, "DELLY2") %>% mutate(SVCaller = "Delly")
svaba_wgscircles = get_svcaller_tx_precision(tx_wgscircles, "SVABA") %>% mutate(SVCaller = "Svaba")
smufin_wgscircles = get_svcaller_tx_precision(tx_wgscircles, "SMUFIN") %>% mutate(SVCaller = "Smufin")
novobreak_wgscircles = get_svcaller_tx_precision(tx_wgscircles, "NOVOBREAK") %>% mutate(SVCaller = "Novobreak")
brass_wgscircles = get_svcaller_tx_precision(tx_wgscircles, "BRASS") %>% mutate(SVCaller = "Brass")
union_wgscircles = tx_wgscircles %>% mutate(SVCaller = "Union")
atleasttwo_wgscircles = tx_wgscircles %>% filter((Smufin+Delly+Svaba+Brass+Novobreak)>=2)  %>% mutate(SVCaller = "AtLeastTwo")
tx_wgscircles = bind_rows(
  delly_wgscircles %>% dplyr::select(-Filter, -ID, -Orientation),
  smufin_wgscircles %>% dplyr::select(-Filter, -ID, -Orientation),
  svaba_wgscircles %>% dplyr::select(-Filter, -ID, -Orientation),
  novobreak_wgscircles %>% dplyr::select(-Filter, -ID, -Orientation),
  brass_wgscircles %>% dplyr::select(-Filter, -ID, -Orientation),
  union_wgscircles,
  atleasttwo_wgscircles
)

tx_circles = rbind(tx_circleseq, tx_wgscircles) %>% as_tibble()
tx = tx_circles %>% dplyr::select(Sample, Cohort, SVCaller, ChrA, PosA, ChrB, PosB, 
                                  MergeInfo_Filter, MergeInfo_Filter, MergeInfo_ReadSupport, BPID) %>% distinct()
tx = tx[!duplicated(paste0(tx$SVCaller, ":", tx$BPID)),] # There are duplicates (due to the fact that in the tables like smufin_wgscircles, ... each call has a line but might have the same breakpoitn pair)
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/MergedTx.Rdata")


## Testing July 10

tx %>% .$BPID %>% unique() %>% length() # 1262
tx %>% filter(SVCaller == "Union") %>% .$BPID %>% unique() %>% length() # 1262 -- that is as expected

tx_nonunion = tx %>% filter(SVCaller != "Union")
tx_nonunion[duplicated(tx_nonunion$BPID,fromLast=T) | duplicated(tx_nonunion$BPID,fromLast=F),] %>% View()

tx_union = tx %>% filter(SVCaller == "Union")
tx_union[duplicated(tx_union$BPID,fromLast=T) | duplicated(tx_union$BPID,fromLast=F),] %>% View()

nrow(tx)
tx_union = tx %>% filter(SVCaller == "Union")
tx_union[duplicated(tx_union$BPID,fromLast=T) | duplicated(tx_union$BPID,fromLast=F),] %>% nrow()
nrow(tx)
tx_union = tx %>% filter(SVCaller == "Union")
tx_union[duplicated(tx_union$BPID,fromLast=T) | duplicated(tx_union$BPID,fromLast=F),] %>% nrow()

## PRECISION INFORMATION
# tx_precision = read.table("~/Desktop/PalmTrees/Data/EliasMail190704/ALL_PATIENTS_merge_translocations_diffvariantcallers_noVCFdups_collapsedups_PASS_filtering_bas99_mapq20_hg19_infostrands_precisioninfo_3_7_19.txt",
#                           sep="\t",
#                           header=F) %>% as_tibble()
# colnames(tx_precision) = c("Sample", "ChrA", "PosA", "ChrB", "PosB", "SVCallerInfo", "MergeInfo")
# tx_precision = tx_precision %>%
#   mutate(Sample = gsub("CB", "NB", Sample),
#          Cohort = ifelse(grepl("NBL", Sample), "Peifer", "Berlin"),
#          ChrA = paste0("chr", ChrA),
#          PosA = as.integer(PosA),
#          ChrB = paste0("chr", ChrB), 
#          PosB = as.integer(PosB)) %>%   
#   separate(MergeInfo, into=c("MergeInfo_Filter", "MergeInfo_nbkps", "MergeInfo_ReadSupport"), sep=",") %>%
#   mutate(Svaba = grepl("SVABA", SVCallerInfo),
#          Smufin = grepl("SMUFIN", SVCallerInfo),
#          Delly = grepl("DELLY2", SVCallerInfo),
#          Brass = grepl("BRASS", SVCallerInfo),
#          Novobreak = grepl("NOVOBREAK", SVCallerInfo)) %>% 
#   filter(Sample != "NBL50", Sample != "NBL49", Sample != "NBL47") %>%
#   filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
# tx_precision$BPID = mapply(makeBPID, tx_precision$Sample, tx_precision$ChrA, tx_precision$PosA, tx_precision$ChrB, tx_precision$PosB)
# 
# delly_precision = get_svcaller_tx_precision(tx_precision, "DELLY2") %>% mutate(SVCaller = "Delly") 
# svaba_precision = get_svcaller_tx_precision(tx_precision, "SVABA") %>% mutate(SVCaller = "Svaba")
# smufin_precision = get_svcaller_tx_precision(tx_precision, "SMUFIN") %>% mutate(SVCaller = "Smufin")
# novobreak_precision = get_svcaller_tx_precision(tx_precision, "NOVOBREAK") %>% mutate(SVCaller = "Novobreak")
# brass_precision = get_svcaller_tx_precision(tx_precision, "BRASS") %>% mutate(SVCaller = "Brass")
# union_precision = tx_precision %>% mutate(SVCaller = "Union")
# atleasttwo_precision = tx_precision %>% filter((Smufin+Delly+Svaba+Brass+Novobreak)>=2)  %>% mutate(SVCaller = "AtLeastTwo")
# 

# # How many duplicates are there in the datasets?
# delly_precision = delly_precision %>% mutate(BP = paste0(ChrA, ":", PosA, ">", ChrB, ":", PosB))
# delly_precision[duplicated(delly_precision$BP, fromLast = T) | duplicated(delly_precision$BP, fromLast=F),] %>% nrow()
# smufin_precision = smufin_precision %>% mutate(BP = paste0(ChrA, ":", PosA, ">", ChrB, ":", PosB))
# smufin_precision[duplicated(smufin_precision$BP, fromLast = T) | duplicated(smufin_precision$BP, fromLast=F),] %>% nrow()
# svaba_precision = svaba_precision %>% mutate(BP = paste0(ChrA, ":", PosA, ">", ChrB, ":", PosB))
# svaba_precision[duplicated(svaba_precision$BP, fromLast = T) | duplicated(svaba_precision$BP, fromLast=F),] %>% nrow()
# brass_precision = brass_precision %>% mutate(BP = paste0(ChrA, ":", PosA, ">", ChrB, ":", PosB))
# brass_precision[duplicated(brass_precision$BP, fromLast = T) | duplicated(brass_precision$BP, fromLast=F),] %>% nrow()
# novobreak_precision = novobreak_precision %>% mutate(BP = paste0(ChrA, ":", PosA, ">", ChrB, ":", PosB))
# novobreak_precision[duplicated(novobreak_precision$BP, fromLast = T) | duplicated(novobreak_precision$BP, fromLast=F),] %>% nrow()
# 
# # 
# 
# inCorrectChrOrder = function(ChrA, PosA, ChrB, PosB){
#   chrord = c("1","2","3","4","5","6","7","8","9","10","11",
#              "12","13","14","15","16","17","18","19","20",
#              "21","22","X","Y","M") %>% paste0("chr",.)
#   if (!(ChrA %in% chrord)) stop("ChrA not a chromsome")
#   if (!(ChrB %in% chrord)) stop("ChrB not a chromsome")
#   if (is.na(ChrA) | is.na(PosA) | is.na(ChrB) | is.na(PosB)) stop("NA not allowed as function input")
#   if (ChrA == ChrB){
#     return(PosA <= PosB)
#   } else {
#     return(which(ChrA == chrord) < which(ChrB == chrord))
#   }
# }
# # inCorrectChrOrder("chr3", 12345678, "chr5", 1) == TRUE
# # inCorrectChrOrder("chr5", 12345678, "chr3", 1) == FALSE
# # inCorrectChrOrder("chr5", 12345678, "chr3", 1) == FALSE
# # inCorrectChrOrder("chrX", 12345678, "chr3", 1) == FALSE
# # inCorrectChrOrder("chrX", 12345678, "chrM", 1) == TRUE
# # inCorrectChrOrder("chrX", 12345678, "chrX", 12345678) == TRUE
# # inCorrectChrOrder("chrX", 12345679, "chrX", 12345678) == FALSE
# tx$inCorrectChrOrder = mapply(inCorrectChrOrder, tx$ChrA, tx$PosA, tx$ChrB, tx$PosB)




