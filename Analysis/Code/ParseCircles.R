rm(list=ls())
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

wgs_ecc = read.table("~/Desktop/PalmTrees/Data/Emails190627/ALL_patients_WGS_inferredcircles_eccDNA.txt")
colnames(wgs_ecc) = c("Chr", "Start", "End", "Sample", "CircleClass")
wgs_ecc$Sample = gsub("CB2", "NB2", wgs_ecc$Sample)
wgs_ecc = as_tibble(wgs_ecc)
wgs_ecc$Method = "WGSCircles"
wgs_ecc = wgs_ecc %>% distinct()

wgs_ec = read.table("~/Desktop/PalmTrees/Data/Emails190627/ALL_patients_WGS_inferredcircles_ecDNA.txt")
colnames(wgs_ec) = c("Chr", "Start", "End", "Sample", "CircleClass")
wgs_ec$Sample = gsub("CB2", "NB2", wgs_ec$Sample)
wgs_ec = as_tibble(wgs_ec)
wgs_ec$Method = "WGSCircles"
wgs_ec = wgs_ec %>% distinct()


circleseq = read.table("~/Desktop/PalmTrees/Data/Emails190627/CircleSeq_ecDNA_eccDNA.bed")
colnames(circleseq) = c("Chr", "Start", "End", "ID", "Strand", "Length", "CircleReads", 
                        "ASCAT_Chr", "ASCAT_Start","ASCAT_End", "ASCAT_NormalTotal", "ASCAT_NormalMinor", "ASCAT_TumorTotal", "ASCAT_TumorMinor", "ASCAT_Length", 
                        "Sample", "CircleClass")
circleseq$Method = "CircleSeq"
circleseq = as_tibble(circleseq) 
#sum(duplicated(circleseq$ID)) # 89
# circleseq[duplicated(circleseq$ID,, fromLast=FALSE) | duplicated(circleseq$ID, fromLast=TRUE),] %>% View()
circleseq = circleseq %>% dplyr::select(Chr, Start, End, Sample, CircleClass, Method)
circleseq = circleseq %>% distinct()
circleseq$Sample = gsub("CB2", "NB2", circleseq$Sample)
circles = bind_rows(wgs_ecc, wgs_ec, circleseq)
circles$Chr = ifelse(circles$Chr == "chrMT", "chrM", circles$Chr)
circles_gr = makeGRangesFromDataFrame(circles,
                                      keep.extra.columns = T,
                                      seqnames = "Chr", start.field = "Start", end.field = "End",
                                      seqinfo = keepStandardChromosomes(seqinfo(BSgenome.Hsapiens.UCSC.hg19)))

# sum(is.na(circles$Chr))
# sum(is.na(circles$Start))
# sum(is.na(circles$End))
# sum(is.na(circles$Sample))
# sum(is.na(circles$CircleClass))
# sum(is.na(circles$Method))
#sum(duplicated(circles %>% mutate(x = paste0(Chr, Start, End, Sample, CircleClass, Method)) %>% .$x))

source("~/Desktop/PalmTrees/Analysis/Code/AlmostDuplicateBreakpoints.R")

cscircles_ecc = circles %>%
  filter(Method == "CircleSeq", CircleClass == "eccDNA") %>%
  mutate(CircleChr = Chr, 
         CircleStart = Start,
         CircleEnd = End) %>%
  dplyr::select(-Chr, -Start, -End)

cscircles_ec = circles %>%
  filter(Method == "CircleSeq", CircleClass == "ecDNA") %>% 
  mutate(CircleChr = Chr, 
         CircleStart = Start,
         CircleEnd = End) %>%
  dplyr::select(-Chr, -Start, -End)

cscircles_ecc$id = paste0(cscircles_ecc$Sample,"_",cscircles_ecc$CircleChr,"_",cscircles_ecc$CircleStart,"_",cscircles_ecc$CircleEnd)
cscircles_ec$id = paste0(cscircles_ec$Sample,"_",cscircles_ec$CircleChr,"_",cscircles_ec$CircleStart,"_",cscircles_ec$CircleEnd)
print(intersect(cscircles_ecc$id, cscircles_ec$id))
cscircles_ecc = cscircles_ecc %>% filter(!(id %in% intersect(cscircles_ecc$id, cscircles_ec$id)))

### ecDNA
wgscircles_ec = circles %>%
  filter(Method == "WGSCircles", CircleClass == "ecDNA") %>% 
  mutate(CircleChr = Chr, 
         CircleStart = Start,
         CircleEnd = End) %>%
  dplyr::select(-Chr, -Start, -End)

merged_wgscircles_ec = list()
samples = unique(wgscircles_ec$Sample)
merged_wgscircles_ec = lapply(samples, function (sample){
  this_circles = wgscircles_ec %>% filter(Sample == sample)
  if (nrow(this_circles)==0) return(data.frame(CircleChr = NA, CircleStart = NA, CircleEnd = NA, CircleLength = NA, Sample = sample))
  this_circles = makeGRangesFromDataFrame(df = this_circles,
                                          seqnames.field = "CircleChr",
                                          start.field = "CircleStart",
                                          end.field = "CircleEnd",
                                          keep.extra.columns = T)
  this_circles = data.frame(reduce(this_circles)) %>%
    dplyr::select(seqnames, start, end, width)
  colnames(this_circles) = c("CircleChr", "CircleStart", "CircleEnd", "CircleLength")
  this_circles$Sample = sample
  return(this_circles)
})
merged_wgscircles_ec = do.call(rbind, merged_wgscircles_ec)
merged_wgscircles_ec = merged_wgscircles_ec %>% filter(isdefchrom(CircleChr)) %>% droplevels()

# WGS eccDNA
wgscircles_ecc = circles %>%
  filter(Method == "WGSCircles", CircleClass == "eccDNA") %>% 
  mutate(CircleChr = Chr, 
         CircleStart = Start,
         CircleEnd = End) %>%
  dplyr::select(-Chr, -Start, -End)

wgscircles_ecc$id = paste0(wgscircles_ecc$Sample,"_",wgscircles_ecc$CircleChr,"_",wgscircles_ecc$CircleStart,"_",wgscircles_ecc$CircleEnd)
wgscircles_ec$id = paste0(wgscircles_ec$Sample,"_",wgscircles_ec$CircleChr,"_",wgscircles_ec$CircleStart,"_",wgscircles_ec$CircleEnd)
print(intersect(wgscircles_ecc$id, wgscircles_ec$id))
wgscircles_ecc = wgscircles_ecc %>% filter(!(id %in% intersect(wgscircles_ecc$id, wgscircles_ec$id)))

merged_wgscircles_ecc = list()
samples = unique(wgscircles_ecc$Sample)
merged_wgscircles_ecc = lapply(samples, function (sample){
  this_circles = wgscircles_ecc %>% filter(Sample == sample)
  if (nrow(this_circles)==0) return(data.frame(CircleChr = NA, CircleStart = NA, CircleEnd = NA, CircleLength = NA, Sample = sample))
  this_circles = makeGRangesFromDataFrame(df = this_circles,
                                          seqnames.field = "CircleChr",
                                          start.field = "CircleStart",
                                          end.field = "CircleEnd",
                                          keep.extra.columns = T)
  this_circles = data.frame(reduce(this_circles)) %>%
    dplyr::select(seqnames, start, end, width)
  colnames(this_circles) = c("CircleChr", "CircleStart", "CircleEnd", "CircleLength")
  this_circles$Sample = sample
  return(this_circles)
})
merged_wgscircles_ecc = do.call(rbind, merged_wgscircles_ecc)
merged_wgscircles_ecc = merged_wgscircles_ecc %>% filter(isdefchrom(CircleChr)) %>% droplevels()


### ecDNA
wgscircles_ec = circles %>%
  filter(Method == "WGSCircles", CircleClass == "ecDNA") %>% 
  mutate(CircleChr = Chr, 
         CircleStart = Start,
         CircleEnd = End) %>%
  dplyr::select(-Chr, -Start, -End)

merged_wgscircles_ec = list()
samples = unique(wgscircles_ec$Sample)
merged_wgscircles_ec = lapply(samples, function (sample){
  this_circles = wgscircles_ec %>% filter(Sample == sample)
  if (nrow(this_circles)==0) return(data.frame(CircleChr = NA, CircleStart = NA, CircleEnd = NA, CircleLength = NA, Sample = sample))
  this_circles = makeGRangesFromDataFrame(df = this_circles,
                                          seqnames.field = "CircleChr",
                                          start.field = "CircleStart",
                                          end.field = "CircleEnd",
                                          keep.extra.columns = T)
  this_circles = data.frame(reduce(this_circles)) %>%
    dplyr::select(seqnames, start, end, width)
  colnames(this_circles) = c("CircleChr", "CircleStart", "CircleEnd", "CircleLength")
  this_circles$Sample = sample
  return(this_circles)
})
merged_wgscircles_ec = do.call(rbind, merged_wgscircles_ec)
merged_wgscircles_ec = merged_wgscircles_ec %>% filter(isdefchrom(CircleChr)) %>% droplevels()

rm(list=c("circleseq", "wgs_ec", "wgs_ecc"))
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/Circles.Rdata")


