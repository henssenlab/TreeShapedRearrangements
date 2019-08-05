rm(list=ls())
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Circles.Rdata") 
library(GenomicRanges)
library(dplyr)
library(tidyr)
source("~/Desktop/PalmTrees/Analysis/Code/AlmostDuplicateBreakpoints.R")

cscircles_ecc = circles %>%
  filter(Method == "CircleSeq", CircleClass == "eccDNA")
cscircles_ec = circles %>%
  filter(Method == "CircleSeq", CircleClass == "ecDNA")

wgscircles_ecc = circles %>%
  filter(Method == "WGSCircles", CircleClass == "eccDNA")
merged_wgscircles_ecc = list()
samples = unique(wgscircles_ecc$Sample)
merged_wgscircles_ecc = lapply(samples, function (sample){
  this_circles = wgscircles_ecc %>% filter(Sample == sample)
  if (nrow(this_circles)==0) return(data.frame(MergedCirclesChr = NA, MergedCirclesStart = NA, MergedCirclesEnd = NA, MergedCirclesLength = NA, Sample = sample))
  this_circles = makeGRangesFromDataFrame(df = this_circles,
                                          seqnames.field = "Chr",
                                          start.field = "Start",
                                          end.field = "End",
                                          keep.extra.columns = T)
  this_circles = data.frame(reduce(this_circles)) %>%
    dplyr::select(seqnames, start, end, width)
  colnames(this_circles) = c("MergedCirclesChr", "MergedCirclesStart", "MergedCirclesEnd", "MergedCirclesLength")
  this_circles$Sample = sample
  return(this_circles)
})
merged_wgscircles_ecc = do.call(rbind, merged_wgscircles_ecc)
merged_wgscircles_ecc = merged_wgscircles_ecc %>% filter(isdefchrom(MergedCirclesChr)) %>% droplevels()


### ecDNA
wgscircles_ec = circles %>%
  filter(Method == "WGSCircles", CircleClass == "ecDNA")
merged_wgscircles_ec = list()
samples = unique(wgscircles_ec$Sample)
merged_wgscircles_ec = lapply(samples, function (sample){
  this_circles = wgscircles_ec %>% filter(Sample == sample)
  if (nrow(this_circles)==0) return(data.frame(MergedCirclesChr = NA, MergedCirclesStart = NA, MergedCirclesEnd = NA, MergedCirclesLength = NA, Sample = sample))
  this_circles = makeGRangesFromDataFrame(df = this_circles,
                                          seqnames.field = "Chr",
                                          start.field = "Start",
                                          end.field = "End",
                                          keep.extra.columns = T)
  this_circles = data.frame(reduce(this_circles)) %>%
    dplyr::select(seqnames, start, end, width)
  colnames(this_circles) = c("MergedCirclesChr", "MergedCirclesStart", "MergedCirclesEnd", "MergedCirclesLength")
  this_circles$Sample = sample
  return(this_circles)
})
merged_wgscircles_ec = do.call(rbind, merged_wgscircles_ec)
merged_wgscircles_ec = merged_wgscircles_ec %>% filter(isdefchrom(MergedCirclesChr)) %>% droplevels()

