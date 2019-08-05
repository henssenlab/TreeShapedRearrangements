library(dplyr)
# source("~/Desktop/PalmTrees/Analysis/Code/ParseSVCallerData.R")
# 
# #system("gzip -d ~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/*/*/*.vcf.gz")
# 
# datapath = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/"
# samples = dir(datapath)
# samples = samples[!grepl("excluded", samples)]
# 
# for (sample in samples){
#   delly = parse_delly_allsv(paste0(datapath, sample, "/delly2/delly2-svs.pre.vcf"), sample)
#   
#   smufin_snvs = parse_smufin(paste0(datapath, sample, "/smufin/somatic_SNVs.txt.sample_", sample), sample)
#   smufin_small = parse_smufin(paste0(datapath, sample, "/smufin/somatic_small_SVs.sample_", sample), sample)
#   smufin_large = parse_smufin(paste0(datapath, sample, "/smufin/somatic_large_SVs.sample_", sample), sample)
# 
#   # Svaba
#   paste0(datapath, sample, "/svaba/svaba-indels.PASS.vcf")
#   paste0(datapath, sample, "/svaba/svaba-sv-fromindels.PASS.vcf")
#   paste0(datapath, sample, "/svaba/svaba-sv.PASS.vcf")
#   
#   # Brass
#   paste0(datapath, sample, "/brass/T_vs_N.annot.vcf")
#   
# }
# 

parse_delly_allsv = function(fname, samplename){
  #fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/CB2003/delly2/delly2-svs.pre.vcf"
  #samplename = "CB2003"
  if (!file.exists(fname)){
    #return(as_tibble(data.frame(Sample = samplename, SVCaller= "Delly", Type= NA, ChrA= NA, PosA= NA, ChrB= NA, PosB= NA, Quality= NA, Precision= NA, PairedEndSupport= NA, PairedEndMedianMAPQ=NA, SplitReadSupport= NA, ConnectionType= NA, DirectionA=NA, DirectionB = NA, JunctionType = NA, SomaticOrGermline= NA)))
    return(NULL)
  }
  v = read.table(fname, header=F, sep="\t", comment.char = "#") %>% as_tibble()
  v$ChrA = paste0("chr", as.character(v$V1))
  v$PosA = v$V2
  v$Type = v$V8 %>% as.character() %>% gsub("^.*SVTYPE=", "", .) %>% gsub(";.*", "", .)
  v$Quality = v$V7 %>% as.character()
  v$Precision = v$V8 %>% as.character() %>% gsub(";.*", "", .)
  v$ChrB = v$V8 %>% gsub("^.*;CHR2=","",.) %>% gsub(";.*", "", .) %>% as.character() %>% paste0("chr", .)
  v$PosB = v$V8 %>% gsub("^.*;END=","",.) %>% gsub(";.*", "", .) %>% as.numeric()
  v$PairedEndSupport = v$V8 %>% gsub("^.*;PE=","",.) %>% gsub(";.*", "", .) %>% as.numeric() # Paired-end support of the structural variant
  v$PairedEndMedianMAPQ = v$V8 %>% gsub("^.*;MAPQ=","",.) %>% gsub(";.*", "", .) %>% as.numeric() # Median mapping quality of paired-ends
  v$SplitReadSupport = v$V8 %>% gsub("^.*;SR=","",.) %>% gsub(";.*", "", .) %>% as.numeric() # Median mapping quality of paired-ends
  v$SplitReadSupport = ifelse(is.na(v$SplitReadSupport), 0, v$SplitReadSupport)
  v$ConnectionType = v$V8 %>% gsub("^.*;CT=","",.) %>% gsub(";.*", "", .) 
  
  v$DirectionA = ifelse(grepl("3to", v$ConnectionType), "Tail", ifelse(grepl("5to", v$ConnectionType), "Head", NA))
  v$DirectionB = ifelse(grepl("to3", v$ConnectionType), "Tail", ifelse(grepl("to5", v$ConnectionType), "Head", NA))
  v$JunctionType = paste0(v$DirectionA, ">", v$DirectionB)
  
  v$SomaticOrGermline = v$V8 %>% gsub("^.*;","",.)
  v$Sample = samplename %>% gsub("CB","NB", .)
  v$SVCaller = "Delly"
  v = v %>% dplyr::select(Sample, SVCaller, Type, ChrA, PosA, ChrB, PosB, Quality, Precision, PairedEndSupport, PairedEndMedianMAPQ, SplitReadSupport, ConnectionType, DirectionA, DirectionB, JunctionType, SomaticOrGermline)
  return(v)
}

parse_smufin_allsv = function(fname_smallSV, fname_largeSV, samplename){
  #fname_smallSV = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/CB2004/smufin/somatic_small_SVs.sample_CB2004"
  #fname_largeSV = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/CB2004/smufin/somatic_large_SVs.sample_CB2004"
  #samplename = "CB2004"
  small = read.table(fname_smallSV, sep=" ", skip=1) %>% as_tibble()
  small$Type = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .) %>% gsub("\t.*", "", .)
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t[[:alnum:]]*\t", "", .)
  small$ChrA = small$V1 %>% gsub("\t.*", "", .) %>% paste0("chr", .)
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$PosA = small$V1 %>% gsub("\t.*", "", .) %>% as.numeric()
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$Length = small$V1 %>% gsub("\t.*", "", .) %>% as.numeric()
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$InsertedSequence = ifelse(small$Type == "INS", small$V1 %>% gsub("\t.*", "", .), NA)
  small$V1 = ifelse(small$Type == "INS", small$V1 %>% gsub("^[[:alnum:]]*\t", "", .), small$V1)
  small$Smufin.NinControl = small$V1 %>% gsub("\t.*", "", .) %>% as.numeric()
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$Smufin.TinControl = small$V1 %>% gsub("\t.*", "", .) %>% as.numeric()
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$Smufin.NinTumor = small$V1 %>% gsub("\t.*", "", .) %>% as.numeric()
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$Smufin.TinTumor = small$V1 %>% gsub("\t.*", "", .) %>% as.numeric()
  small$V1 = small$V1 %>% gsub("^[[:alnum:]]*\t", "", .)
  small$ChrB = small$ChrA
  small$PosB = ifelse((small$Type == "INS"), (small$PosA + 1), (small$PosA + small$Length + 1))
  small = small %>% dplyr::select(-V1)
  
  large = read.table(fname_largeSV, header=T, sep="\t") %>% as_tibble()
  large$Type = ifelse(large$Type == "BKP", "BND", NA)
  large$ChrA = paste0("chr", large$Chr_BKP_1)
  large$PosA = large$Pos_BKP_1 %>% as.numeric()
  large$ChrB = paste0("chr", large$Chr_BKP_2)
  large$PosB = large$Pos_BKP_2 %>% as.numeric()
  large$Smufin.BNDCoverage = large$Coverage %>% as.numeric()
  large = large %>% dplyr::select(Type, ChrA, PosA, ChrB, PosB, Smufin.BNDCoverage)
  
  small$Smufin.BNDCoverage = NA
  large$Length = NA
  large$InsertedSequence = NA
  large$Smufin.NinControl = NA
  large$Smufin.TinControl = NA
  large$Smufin.NinTumor = NA
  large$Smufin.TinTumor = NA
  
  all = bind_rows(small, large)
  all$Sample = samplename %>% gsub("CB","NB", .)
  all$SVCaller = "Smufin"
  all = all %>% dplyr::select(Sample, SVCaller, Type, ChrA, PosA, ChrB, PosB, 
                              Smufin.NinControl, Smufin.TinControl, Smufin.NinTumor, Smufin.TinTumor,
                              Smufin.BNDCoverage)
}

parse_svaba_allsv = function(indels_fname, svfromindels_fname, sv_fname, sample=NA, cohort=NA){
  
  # indels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2013/svaba/svaba-indels.PASS.vcf"
  # svfromindels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2013/svaba/svaba-sv-fromindels.PASS.vcf"
  # sv_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2013/svaba/svaba-sv.PASS.vcf"
  # sample = "NB2013"
  # cohort=NA

  # indels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2001/svaba/svaba-indels.PASS.vcf"
  # svfromindels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2001/svaba/svaba-sv-fromindels.PASS.vcf"
  # sv_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2001/svaba/svaba-sv.PASS.vcf"
  # sample = "NB2001"
  # cohort=NA
  
  # indels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/CB2027/svaba/svaba-indels.PASS.vcf"
  # svfromindels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/CB2027/svaba/svaba-sv-fromindels.PASS.vcf"
  # sv_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/CB2027/svaba/svaba-sv.PASS.vcf"
  # sample = "CB2027"
  # cohort=NA
  
  indels_vcf = tryCatch(read.table(indels_fname,
                                   comment.char = "#",
                                   header = F,
                                   sep = "\t") %>% as_tibble(),
                        error = function(e) NULL)
  
  if (!is.null(indels_vcf)){
    colnames(indels_vcf) = c("ChrA", "PosA", "ID","REF","ALT","QUAL","FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR")
    indels_vcf$Sample = sample
    indels_vcf$Cohort = cohort
    indels_vcf$Class = "indel"
    indels_vcf$ChrA = as.character(indels_vcf$ChrA)
    indels_vcf$PosA = as.numeric(indels_vcf$PosA)
    indels_vcf$ChrB = indels_vcf$ChrA
    indels_vcf$PosB = indels_vcf$PosA + nchar(as.character(indels_vcf$REF)) - 1 # checked, those are exactly the ucsc
    indels_vcf$Filter = indels_vcf$FILTER
    
    indels_vcf = indels_vcf %>% 
      dplyr::select(Sample, Cohort, Class, ChrA, PosA, ChrB, PosB, Filter)
  }
  
  ############################
  
  svfromindels_vcf = tryCatch(read.table(svfromindels_fname,
                                         comment.char = "#",
                                         header = F,
                                         sep = "\t") %>% as_tibble(),
                              error = function(e) NULL)
  
  if (!is.null(svfromindels_vcf)){
    colnames(svfromindels_vcf) = c("ChrA",
                                   "PosA",
                                   "ID",
                                   "REF",
                                   "ALT",
                                   "QUAL",
                                   "FILTER",
                                   "INFO",
                                   "FORMAT",
                                   "NORMAL", 
                                   "TUMOR")
    svfromindels_vcf$Sample = sample
    svfromindels_vcf$Cohort = cohort
    svfromindels_vcf$Class = "indelfromsv"
    svfromindels_vcf$ChrA = as.character(svfromindels_vcf$ChrA)
    svfromindels_vcf$PosA = as.numeric(svfromindels_vcf$PosA)
    svfromindels_vcf$ChrB = svfromindels_vcf$ChrA
    svfromindels_vcf$PosB = svfromindels_vcf$PosA + nchar(as.character(svfromindels_vcf$REF)) - 1
    svfromindels_vcf = svfromindels_vcf %>% 
      mutate(Filter=FILTER) %>%
      dplyr::select(Sample, Cohort, Class, ChrA, PosA, ChrB, PosB, Filter)
  }
  
  ############################
  
  sv_vcf = tryCatch(read.table(sv_fname,
                               comment.char = "#",
                               header = F,
                               sep = "\t") %>% as_tibble(),
                    error = function(e) NULL)

  if (!is.null(sv_vcf)){
    
    colnames(sv_vcf) = c("ChrA","PosA","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL", "TUMOR")
    
    sv_vcf$BNDPairID = sv_vcf$ID %>% gsub(":.*$", "", .)
    sv_vcf$ChrA = as.character(sv_vcf$ChrA)
    sv_vcf$PosA = as.numeric(sv_vcf$PosA)
    sv_vcf$Sample = sample
    sv_vcf$Type = sapply(sv_vcf$ID, function(ID) svaba_get_sv_type(ID, sv_vcf)) 
    sv_vcf$BracketCoord = sv_vcf$ALT
    sv_vcf$DirectionA = ifelse(grepl("^[AGCTN]", sv_vcf$ALT), "Tail", "Head") # Head = Breakpoint Coordinate is Max Coordinate = Coordinate decrease with distance from Breakpoint
    sv_vcf$DirectionB = ifelse(grepl("\\[", sv_vcf$ALT), "Head", "Tail")
    sv_vcf$JunctionType = paste0(sv_vcf$DirectionA, ">", sv_vcf$DirectionB)

    sv_vcf$ALT = ifelse(
      grepl("\\[", sv_vcf$ALT),
      sv_vcf$ALT %>% gsub("\\[[[:alnum:]]*$", "", .) %>% gsub("^[[:alnum:]]*\\[", "", .),
      sv_vcf$ALT %>% gsub("\\][[:alnum:]]*$", "", .) %>% gsub("^[[:alnum:]]*\\]", "", .)
    )
    
    sv_vcf$ChrB = gsub(":.*$", "", sv_vcf$ALT) %>% as.character() # some go wrong, because
    sv_vcf$PosB = gsub("^.*:", "", sv_vcf$ALT) %>% as.numeric() # some go wrong,
    sv_vcf$Evidence = sv_vcf$INFO %>% gsub("^.*EVDNC=", "", .) %>% gsub(";.*", "", .) %>% as.character()
    sv_vcf$Homology = ifelse(
      grepl("HOMSEQ", sv_vcf$INFO),
      sv_vcf$INFO %>% gsub("^.*;HOMSEQ=", "", .) %>% gsub(";.*", "", .) %>% as.character(),
      ""
    )
    sv_vcf$Homology = ifelse(grepl(";IMPRECISE", sv_vcf$INFO),
                             NA,
                             sv_vcf$Homology)
    sv_vcf$HomologyLength = nchar(sv_vcf$Homology)

    sv_vcf$Insertion = ifelse(
      grepl("INSERTION", sv_vcf$INFO),
      sv_vcf$INFO %>% gsub("^.*;INSERTION=", "", .) %>% gsub(";.*", "", .) %>% as.character(),
      ""
    )
    sv_vcf$Insertion = ifelse(grepl(";IMPRECISE", sv_vcf$INFO),
                             NA,
                             sv_vcf$Insertion)
    sv_vcf$InsertionLength = nchar(sv_vcf$Insertion)

    sv_vcf$Class = "sv"
    sv_vcf$Cohort = cohort
    #sv_vcf$Type = sv_vcf$INFO %>% gsub("^.*;SVTYPE=", "", .) %>% gsub(";.*", "", .) %>% as.character()
    sv_vcf$MappingQuality = sv_vcf$INFO %>% gsub("^.*;MAPQ=", "", .) %>% gsub(";.*", "", .) %>% as.numeric()
    sv_vcf$Quality = sv_vcf$QUAL %>% as.numeric()
    
    if (length(unique(sv_vcf$FORMAT)) > 1){
      stop("FORMAT row is not identical for different SVs.")
    } else {
      sv_vcf_format = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$FORMAT), "\\:", perl=T)))
      sv_vcf_normal = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$NORMAL), "\\:", perl=T)))
      colnames(sv_vcf_normal) = paste0("Normal_",sv_vcf_format[1,])
      sv_vcf_tumor = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$TUMOR), "\\:", perl=T)))
      colnames(sv_vcf_tumor) = paste0("Tumor_",sv_vcf_format[1,])
    }
    sv_vcf = cbind(sv_vcf, sv_vcf_normal, sv_vcf_tumor)
    
    sv_vcf = 
      sv_vcf %>%
      filter(grepl(":1$", ID)) # assuming that *:1 for intrachromosomal rearrangements always is PosA < PosB
    
    sv_vcf = sv_vcf %>% mutate(Filter=FILTER, Tumor_LocalCoverage=as.numeric(Tumor_DP), TumorSplitReads=as.numeric(Tumor_SR), TumorDiscordantReads=as.numeric(Tumor_DR), TumorReadSupport=as.numeric(Tumor_AD),
                               NormalSplitReads=as.numeric(Normal_SR), NormalDiscordantReads=as.numeric(Normal_DR), NormalReadSupport=as.numeric(Normal_AD)) %>%
      dplyr::select(Sample, Cohort, Class, BNDPairID, ID, ChrA, PosA, ChrB, PosB, Filter, Type, DirectionA, DirectionB, JunctionType, Evidence, Tumor_LocalCoverage, TumorReadSupport, TumorDiscordantReads, TumorSplitReads, NormalReadSupport, NormalDiscordantReads, NormalSplitReads, Homology,HomologyLength,Insertion,InsertionLength)
  }
  
  vcf = bind_rows(sv_vcf, indels_vcf, svfromindels_vcf) %>% as_tibble()
  
  vcf$ChrA  = vcf$ChrA %>% gsub("chr", "", .) %>% paste0("chr", .)
  vcf$ChrB  = vcf$ChrB %>% gsub("chr", "", .) %>% paste0("chr", .)

  return(vcf)
}

svaba_get_sv_type <- function(svaba_vcf_table_ID, svaba_vcf_table){
  # adapted from: https://github.com/walaj/svaba/issues/4
  
  #svaba_vcf_table_ID = sv_vcf[[1, "ID"]]
  #svaba_vcf_table = sv_vcf
  
  root <- gsub(":[12]", "", svaba_vcf_table_ID)
  mate1 <- paste0(root, ":1")
  mate2 <- paste0(root, ":2")
  
  alt1 = tryCatch(svaba_vcf_table %>% filter(ID == mate1) %>% .$ALT, error = function(e) NA) %>% as.character()
  alt2 = tryCatch(svaba_vcf_table %>% filter(ID == mate2) %>% .$ALT, error = function(e) NA) %>% as.character()
  chr1 = tryCatch(svaba_vcf_table %>% filter(ID == mate1) %>% .$ChrA, error = function(e) NA) %>% as.character()
  chr2 = tryCatch(svaba_vcf_table %>% filter(ID == mate2) %>% .$ChrA, error = function(e) NA) %>% as.character()
  
  if (identical(alt1, character(0)) | identical(alt2, character(0)) | identical(chr1, character(0)) | identical(chr2, character(0))) return(NA)
  if (is.na(alt1) | is.na(alt2) | is.na(chr1) | is.na(chr2)) return(NA)
  
  if ((grepl("\\[", alt1) & grepl("\\[", alt2)) | (grepl("\\]", alt1) & grepl("\\]", alt2))){
    sv_type <- "INV"
    
  } else if (grepl("[ACTGN]\\[", alt1) & grepl("^\\]", alt2)){
    sv_type <- "DEL"
    
  } else if (grepl("^\\]", alt1) & grepl("[ACTGN]\\[", alt2)){
    sv_type <- "DUP/INS"
    
  } else{
    sv_type <- "UNKNOWN"
  }
  
  # own category BND for interchromosomal rearrangements
  if (chr1 != chr2) sv_type = "BND"
  
  return(sv_type)
}

svaba_get_large_insertions = function(sv){
  # Assumes that each breakpoint appears twice in sv table (as in raw svaba output)
  
  # indels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2013/svaba/svaba-indels.PASS.vcf"
  # svfromindels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2013/svaba/svaba-sv-fromindels.PASS.vcf"
  # sv_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NB2013/svaba/svaba-sv.PASS.vcf"
  # sample = "NB2013"
  # cohort=NA
  # 
  # indels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NBL04/svaba/svaba-indels.PASS.vcf"
  # svfromindels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NBL04/svaba/svaba-sv-fromindels.PASS.vcf"
  # sv_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NBL04/svaba/svaba-sv.PASS.vcf"
  # sample = "NBL04"
  # cohort=NA
  
  
  # indels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NBL42/svaba/svaba-indels.PASS.vcf"
  # svfromindels_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NBL42/svaba/svaba-sv-fromindels.PASS.vcf"
  # sv_fname = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/NBL42/svaba/svaba-sv.PASS.vcf"
  # sample = "NBL42"
  # cohort=NA

  #sv = parse_svaba_allsv(indels_fname, svfromindels_fname, sv_fname, sample, cohort)
    
  sv = sv %>% filter(Class == "sv")
  
  insertions = list()
  k=1
  for (i in 1:nrow(sv)){
    
    if (sv[[i, "DirectionA"]] != "Head") next
    
    # ChrA = "host chromosome"; ChrB = "donor chromosome"
    potential_partner_sv = sv %>% filter(BNDPairID != sv[[i,"BNDPairID"]],
                                         ChrA == sv[[i, "ChrA"]], # "A" is the host chromosome where something from ChrB is integrated
                                         sv[[i, "PosA"]] < PosA,
                                         DirectionA == "Tail",
                                         abs(PosA - sv[[i, "PosA"]]) < 1000000, # there is less than ... missing at the potential integration site
                                         ChrB == sv[[i, "ChrB"]], # ChrB is the donor chromosome
                                         abs(PosB - sv[[i, "PosB"]]) < 10000000, # the insertion is smaller than ...
                                         DirectionA != sv[[i,"ChrA"]]) # make sure the directions are good

    if (nrow(potential_partner_sv) == 0) next
    
    potential_partner_sv = 
      potential_partner_sv %>%
      mutate(
        EvidenceLeft = sv[[i,"Evidence"]],
        EvidenceRight = Evidence, 
        ReadSupportLeft = sv[[i,"TumorReadSupport"]],
        ReadSupportRight = TumorReadSupport,
        HostChr = sv[[i, "ChrA"]],
        HostLeft = sv[[i, "PosA"]],
        HostRight = PosA,
        DonorChr = sv[[i, "ChrB"]],
        DonorLeft = sv[[i, "PosB"]],
        DonorRight = PosB
      ) %>%
      dplyr::select(Sample, Cohort, HostChr, HostLeft, HostRight, DonorChr, DonorLeft, DonorRight, ReadSupportLeft, ReadSupportRight) %>% distinct()
    insertions[[k]] = potential_partner_sv
    k = k + 1
  }
  insertions = do.call(rbind, insertions)
  if (!is.null(insertions)) insertions = insertions %>% distinct()
}
