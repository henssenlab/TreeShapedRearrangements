library(dplyr)

# makeBPID = function(Sample, ChrA, PosA, ChrB, PosB){
#   chrord = c("1","2","3","4","5","6","7","8","9","10","11",
#              "12","13","14","15","16","17","18","19","20",
#              "21","22","X","Y","M") %>% paste0("chr",.)
#   if (!(ChrA %in% chrord)) stop("ChrA not a chromsome")
#   if (!(ChrB %in% chrord)) stop("ChrB not a chromsome")
#   if (is.na(ChrA) | is.na(PosA) | is.na(ChrB) | is.na(PosB)) stop("NA not allowed as function input")
#   if (ChrA == ChrB){
#     return(ifelse(PosA <= PosB, paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB), paste0(Sample, "_", ChrB, ":", PosB, "_", ChrA, ":", PosA)))
#   } else {
#     return(ifelse(which(ChrA == chrord) < which(ChrB == chrord), paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB), paste0(Sample, "_", ChrB, ":", PosB, "_", ChrA, ":", PosA)))
#   }
# }

parse_svaba_circlejunction_vcf = function(sv_fname, sample, cohort=NA){

  # sv_fname = "~/Desktop/PalmTrees/Data/CircleSeq_CircleReads_Svaba/AH_01_KELLY_1_val.hg19.sorted.RmDup.CircPairTestF49.QFilt20.merged.bam.svaba.sv.vcf"
  # sample = NA
  # cohort=NA
  # 
  # sv_fname = "~/Desktop/PalmTrees/Data/CircleSeq_CircleReads_Svaba/AH_03_SK-N-SH_1_val.hg19.sorted.RmDup.CircPairTestF49.QFilt20.merged.bam.svaba.sv.vcf"
  # sample = NA
  # cohort=NA

  sv_vcf = tryCatch(read.table(sv_fname,
                               comment.char = "#",
                               header = F,
                               sep = "\t") %>% as_tibble(),
                    error = function(e) NULL)
  
  if (!is.null(sv_vcf)){
    
    colnames(sv_vcf) = c("ChrA","PosA","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","NORMAL")
    
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
    
    sv_vcf$isPrecise = !grepl("IMPRECISE", sv_vcf$INFO) 
    sv_vcf = sv_vcf %>% filter(ChrA == ChrB, isPrecise) # <------------------------------------ We only take intrachromosomal and precise SV calls into account
    
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
    sv_vcf$MappingQuality = sv_vcf$INFO %>% gsub("^.*;MAPQ=", "", .) %>% gsub(";.*", "", .) %>% as.numeric()
    sv_vcf$Quality = sv_vcf$QUAL %>% as.numeric()
    
    if (length(unique(sv_vcf$FORMAT)) > 1){
      stop("FORMAT row is not identical for different SVs.")
    } else {
      sv_vcf_format = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$FORMAT), "\\:", perl=T)))
      sv_vcf_normal = as_tibble(do.call(rbind, strsplit(as.character(sv_vcf$NORMAL), "\\:", perl=T)))
      colnames(sv_vcf_normal) = paste0("Normal_",sv_vcf_format[1,])
    }
    sv_vcf = cbind(sv_vcf, sv_vcf_normal)
    
    sv_vcf = 
      sv_vcf %>%
      filter(grepl(":1$", ID)) # assuming that *:1 for intrachromosomal rearrangements always is PosA < PosB
    
    sv_vcf = sv_vcf %>% mutate(Filter=FILTER, 
                               NormalSplitReads=as.numeric(Normal_SR), NormalDiscordantReads=as.numeric(Normal_DR), NormalReadSupport=as.numeric(Normal_AD)) %>%
      dplyr::select(Sample, Cohort, Class, BNDPairID, ID, ChrA, PosA, ChrB, PosB, Filter, Type, DirectionA, DirectionB, JunctionType, Evidence, NormalReadSupport, NormalDiscordantReads, NormalSplitReads, Homology,HomologyLength,Insertion,InsertionLength,isPrecise)
  }
  
  svaba_circles = sv_vcf %>% filter(Type == "DUP/INS") %>% as_tibble()
  
  svaba_circles$ChrA  = svaba_circles$ChrA %>% gsub("chr", "", .) %>% paste0("chr", .)
  svaba_circles$ChrB  = svaba_circles$ChrB %>% gsub("chr", "", .) %>% paste0("chr", .)
  
  return(svaba_circles)
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
