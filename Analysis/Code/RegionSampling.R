library(tidyr)
library(dplyr)

sample_random_regions = function(n, segment_length){
  hg19idx = read.csv("~/Desktop/eccDNADataAnalysis/hg19test.txt", header=F, sep="\t")
  hg19idx = hg19idx[1:24,]
  hg19idx$p = hg19idx$V2 / sum(hg19idx$V2)
  hg19idx = hg19idx[order(hg19idx$p),]
  # draw chromosome
  regions = data.frame(Chr = base::sample(hg19idx$V1, n, replace=TRUE, prob=hg19idx$p))
  regions$max = as.numeric(lapply(regions$Chr, function(c) hg19idx[hg19idx$V1 == c, "V2"]))  - segment_length
  # draw region
  regions$Start = ceiling(runif(n)*regions$max)
  regions$End = regions$Start + segment_length
  return(regions)
}

sample_random_regions_onechr = function(n, segment_length, chr){
  hg19idx = read.csv("~/Desktop/eccDNADataAnalysis/hg19test.txt", header=F, sep="\t")
  max_start_pos = hg19idx[hg19idx$V1 == chr, "V2"] - segment_length + 1
  regions = data.frame(ID=c(paste0("RandomRegion", as.character(1:n))), Chr=chr, Start=ceiling(runif(n)*max_start_pos))
  regions$End = regions$Start + segment_length - 1
  return(regions)
}

sample_random_regions_broad = function(n, segment_length){
  # randomly samples "n" intervals of length "segment_length" from a masked genome defined in "broad"
  # returns a nx3 data.frame with columns Chr, Start and End
  broad = read.table("~/Desktop/PalmTrees/Data/BroadIntervals/b37-wgs_calling_regions.v1.interval_list_toHg19_chr_MinusBlckLst.txt", sep="\t", header=F)
  broad$ID = paste0("ID_", as.character(1:nrow(broad)))
  broad$Length = broad$V3-broad$V2 # Length of all whitelisted intervals
  broad = broad[(broad$Length > segment_length),] 
  regions = data.frame(ID = base::sample(broad$ID, n, replace=TRUE, prob=broad$Length)) # sample a whitelisted interval, probability according to its length
  regions$max = broad[match(regions$ID, broad$ID), "Length"] - segment_length # maximal start base for the random sample (sample must not exceed whitelist interval)
  regions$Chr = broad[match(regions$ID, broad$ID), "V1"] # chromosome 
  regions$Start = ceiling(runif(n)*regions$max) + as.numeric(broad[match(regions$ID, broad$ID), "V2"]) # sample a start position
  regions$End = regions$Start + segment_length # end position
  regions = regions[,-c(1,2)]
}
