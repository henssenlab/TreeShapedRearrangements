library(tidyr)
library(dplyr)

# sampling method used for all randomization analyses in Koche et al. 
# randomly samples "n" intervals of length "segment_length" from a 
# genome whitelist defined in "broad"
# returns a nx3 data.frame with columns Chr, Start and End
sample_random_regions_broad = function(n, segment_length){
  broad = read.table("~/Desktop/PalmTrees/Data/BroadIntervals/b37-wgs_calling_regions.v1.interval_list_toHg19_chr_MinusBlckLst.txt", sep="\t", header=F)
  broad$ID = paste0("ID_", as.character(1:nrow(broad)))
  broad$Length = broad$V3-broad$V2 # length of all whitelisted intervals
  broad = broad[(broad$Length > segment_length),] 
  regions = data.frame(ID = base::sample(broad$ID, n, replace=TRUE, prob=broad$Length)) # sample a whitelisted interval, probability according to its length
  regions$max = broad[match(regions$ID, broad$ID), "Length"] - segment_length # maximal start base for the random sample (sample must not exceed whitelist interval)
  regions$Chr = broad[match(regions$ID, broad$ID), "V1"] # chromosome 
  regions$Start = ceiling(runif(n)*regions$max) + as.numeric(broad[match(regions$ID, broad$ID), "V2"]) # sample a start position
  regions$End = regions$Start + segment_length # end position
  regions = regions[,-c(1,2)]
}
