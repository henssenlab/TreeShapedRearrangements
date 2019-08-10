rm(list=ls())
library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(Biostrings)
library(regioneR)

broad = read.table("~/Desktop/PalmTrees/Data/BroadIntervals/b37-wgs_calling_regions.v1.interval_list_toHg19_chr_MinusBlckLst.txt", sep="\t", header=F) %>% as_tibble()
broad$V4 = "*" # there are only + in this
colnames(broad) = c("chr", "start", "end", "strand", "interval_name")
broad_gr = makeGRangesFromDataFrame(broad)
broad_mask = gaps(broad_gr, start=1, end=seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)])
broad_mask = broad_mask[strand(broad_mask) == "*"]
broad_mask = c(broad_mask, GRanges("chrY:1-59373566"))
human.genome <- getGenomeAndMask("hg19", mask=NA)$genome
human.canonical = filterChromosomes(human.genome, organism="hg", chr.type="canonical")
# example usage: r=randomizeRegions(this_snv_gr, genome=human.canonical, mask=broad_mask, allow.overlaps=FALSE)
# example usage: r=circularRandomizeRegions(this_snv_gr, genome=human.canonical, maks=broad_mask)

nregions = 1000000
mean_seq_lengths = 61
sd_seq_lengths = 0
BreakpointIntervals.gr = createRandomRegions(nregions = nregions, length.mean=mean_seq_lengths, length.sd = sd_seq_lengths, genome=human.canonical, mask=broad_mask, non.overlapping=F)
seqlevels(BreakpointIntervals.gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(BreakpointIntervals.gr) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
BreakpointIntervals.gr <- trim(BreakpointIntervals.gr)
BreakpointIntervals.gr$ID = paste0("RandomSeqs_", seqnames(BreakpointIntervals.gr), "_", as.character(start(BreakpointIntervals.gr)+floor(seqlengths(BreakpointIntervals.gr)/2)))
seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, BreakpointIntervals.gr, as.character=F)
seq = seq[width(seq)>20]
names(seq) = BreakpointIntervals.gr$ID
writeXStringSet(seq, "~/Desktop/PalmTrees/Data/1MRandomSeq_61bp.fa")