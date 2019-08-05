library(dplyr)
library(tidyr)
library(GenomicRanges)
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

myGetSeq = function (gr) getSeq(BSgenome.Hsapiens.UCSC.hg19, gr, as.character=T)

# get_this_cutoff = function(i){
#   if (i <= 20){
#     return(0)
#   } else if (i <= 100){
#     return(floor(i/5)-3)
#   } else {
#     return(17)  
#   }
# }

get_this_cutoff = function(i){
  if (i <= 10){
    return(0)
  } else if (i <= 15){
    return(1)
  } else if (i <= 20){
    return(2)
  } else if (i <= 25){
    return(3) 
  } else {
    return(4) 
  }
}

getHomologySvabaBreakpoints = function(ChrA, PosA, DirectionA, ChrB, PosB, DirectionB){
  
  if ((ChrA == ChrB) & (PosB<PosA)){
    ChrTmp = ChrA
    PosTmp = PosA
    DirectionTmp = DirectionA
    ChrA=ChrB
    PosA=PosB
    DirectionA=DirectionB
    ChrB=ChrTmp
    PosB=PosTmp
    DirectionB=DirectionTmp
  }
  
  dpos = 150
  
  if ((ChrA == ChrB) & ((PosB-PosA)<50)) return(NA)
  
  # The following is necessary not to be fooled by very small duplications that
  # would otherwise be called as one huge microhomology
  dpos = ifelse(ChrA == ChrB,
                min(dpos, abs(PosB-PosA)-get_this_cutoff(abs(PosB-PosA))),
                dpos)
  dpos = max(0,dpos) ## That prevents negative dpos for very close breakpoints  
  
  JunctionType = paste0(DirectionA, ">", DirectionB)
  this_adist = Inf
  
  `%notin%` <- Negate(`%in%`)
  if (ChrA %notin% standardChromosomes(BSgenome.Hsapiens.UCSC.hg19) |
      ChrB %notin% standardChromosomes(BSgenome.Hsapiens.UCSC.hg19) |
      PosA <= 250 | PosB <= 250){
    return(NA)
  }
  
  if (JunctionType == "Head>Tail"){
    RealA = GRanges(paste0(ChrA,":",PosA,"-",PosA+dpos,":","+")) %>% myGetSeq %>% DNAString
    RealB = GRanges(paste0(ChrB, ":", PosB-dpos, "-", PosB, ":", "+")) %>% myGetSeq %>% DNAString
    lRB = nchar(RealB)
    hom=""
    for (i in 0:dpos){ # It must start from zero, otherwise it does not detect 1bp microhomologies
      sRA = substr(RealA, 1, 1+i)
      sRB = substr(RealB, lRB-i, lRB)
      previous_adist = this_adist
      this_adist = adist(sRA,sRB)
      if (this_adist<=get_this_cutoff(i) & this_adist<=previous_adist) hom=sRA
    }
  } else if (JunctionType == "Tail>Head"){
    RealA = GRanges(paste0(ChrA,":",PosA-dpos,"-",PosA,":","+")) %>% myGetSeq %>% DNAString
    RealB = GRanges(paste0(ChrB,":",PosB,"-",PosB+dpos,":","+")) %>% myGetSeq %>% DNAString
    lRA = nchar(RealA)
    hom = ""
    for (i in 0:dpos){
      sRA = substr(RealA, lRA-i, lRA)
      sRB = substr(RealB, 1, 1+i)
      previous_adist = this_adist
      this_adist = adist(sRA,sRB)
      if (this_adist<=get_this_cutoff(i) & this_adist<=previous_adist) hom=sRA
    }
  } else if (JunctionType == "Head>Head"){
    RealA = GRanges(paste0(ChrA,":",PosA,"-",PosA+dpos,":","+")) %>% myGetSeq %>% DNAString
    RealB = GRanges(paste0(ChrB,":",PosB,"-",PosB+dpos,":","+")) %>% myGetSeq %>% DNAString
    hom=""
    for (i in 0:dpos){
      sRA = substr(RealA, 1, 1+i)
      sRB = reverseComplement(substr(RealB, 1, 1+i))
      previous_adist = this_adist
      this_adist = adist(sRA,sRB)
      if (this_adist<=get_this_cutoff(i) & this_adist<=previous_adist) hom=sRA
    }
  } else if (JunctionType == "Tail>Tail"){
    RealA = GRanges(paste0(ChrA,":",PosA-dpos,"-",PosA,":","+")) %>% myGetSeq %>% DNAString
    RealB = GRanges(paste0(ChrB,":",PosB-dpos,"-",PosB,":","+")) %>% myGetSeq %>% DNAString
    lRA = nchar(RealA)
    lRB = nchar(RealB)
    hom=""
    for (i in 0:dpos){
      sRA = substr(RealA, lRA-i, lRA)
      sRB = reverseComplement(substr(RealB, lRB-i, lRB))
      previous_adist = this_adist
      this_adist = adist(sRA,sRB)
      if (this_adist<=get_this_cutoff(i) & this_adist<=previous_adist) hom=sRA
    }
  } else {
    stop("Invalid DirectionA or DirectionB.")
  }
  return(as.character(hom))
}