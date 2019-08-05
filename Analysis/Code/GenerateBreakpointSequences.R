rm(list=ls())

library(dplyr)
library(tidyr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

precise_BPID = tx_circles %>% filter(grepl(",PRECISE", as.character(SVCallerInfo))) %>% .$BPID %>% as.character() %>% unique()

IntervalLength = 60
svcallers = unique(palmtrees$SVCaller)
for (svcaller in svcallers){
  precise_breakpoints = tx_original %>% filter(SVCaller == svcaller, BPID %in% precise_BPID) %>% dplyr::select(SVCaller, Sample, ChrA, PosA, ChrB, PosB, BPID, Method) %>% distinct()
  palmtreetx_precise = txptinfo %>% filter(SVCaller == svcaller, BPID %in% precise_BPID) %>% dplyr::select(SVCaller, Sample, PalmTreeChrom, PalmTreePos, TargetChrom, TargetPos, BPID) %>% .$BPID %>% as.character() %>% unique()
  precise_breakpoints$isPT = unlist(lapply(precise_breakpoints$BPID, function (s) s %in% palmtreetx_precise))
  precise_breakpoints$PalmTreeOrNoPalmTree = ifelse(precise_breakpoints$isPT, "PalmTree", "NoPalmTree")
  
  print(svcaller)
  print("Number of precise palm tree rearrangements")
  print(sum(precise_breakpoints$isPT))
  print("Percent precise palm tree rearrangements (of all precise rearrangements)")
  print(sum(precise_breakpoints$isPT) / length(precise_breakpoints$isPT)*100)
  
  BreakpointIntervals = data.frame(
    ID = c(paste0(precise_breakpoints$Sample, "_", precise_breakpoints$SVCaller, "_", precise_breakpoints$PalmTreeOrNoPalmTree, "_",precise_breakpoints$ChrA, "_", precise_breakpoints$PosA), 
           paste0(precise_breakpoints$Sample, "_", precise_breakpoints$SVCaller, "_", precise_breakpoints$PalmTreeOrNoPalmTree, "_",precise_breakpoints$ChrB, "_", precise_breakpoints$PosB)),
    chr = c(precise_breakpoints$ChrA, precise_breakpoints$ChrB),
    start = c(precise_breakpoints$PosA - ceiling(IntervalLength / 2), precise_breakpoints$PosB - ceiling(IntervalLength / 2)),
    end = c(precise_breakpoints$PosA + ceiling(IntervalLength / 2), precise_breakpoints$PosB + ceiling(IntervalLength / 2)))
  BreakpointIntervals$strand = "*"
  BreakpointIntervals = BreakpointIntervals %>% distinct()
  
  # All Breakpoints for this caller
  BreakpointIntervals.gr = makeGRangesFromDataFrame(df = BreakpointIntervals, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  seqlevels(BreakpointIntervals.gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(BreakpointIntervals.gr) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  BreakpointIntervals.gr <- trim(BreakpointIntervals.gr)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, BreakpointIntervals.gr, as.character=F)
  seq = seq[width(seq)>20]
  names(seq) = BreakpointIntervals.gr$ID
  writeXStringSet(seq, paste0("~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/", svcaller, "_PTvsNoPT_BreakpointSequences.fa"))
  
  # Only Palm Trees
  BreakpointIntervals.gr = makeGRangesFromDataFrame(df = BreakpointIntervals, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  BreakpointIntervals.gr = BreakpointIntervals.gr[!grepl("NoPalmTree", BreakpointIntervals.gr$ID)]
  seqlevels(BreakpointIntervals.gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(BreakpointIntervals.gr) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  BreakpointIntervals.gr <- trim(BreakpointIntervals.gr)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, BreakpointIntervals.gr, as.character=F)
  seq = seq[width(seq)>20]
  names(seq) = BreakpointIntervals.gr$ID
  writeXStringSet(seq, paste0("~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/", svcaller, "_OnlyPalmTreeTx_BreakpointSequences.fa"))
  
  # Only Non Palm Trees
  BreakpointIntervals.gr = makeGRangesFromDataFrame(df = BreakpointIntervals, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  BreakpointIntervals.gr = BreakpointIntervals.gr[grepl("NoPalmTree", BreakpointIntervals.gr$ID)]
  seqlevels(BreakpointIntervals.gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(BreakpointIntervals.gr) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  BreakpointIntervals.gr <- trim(BreakpointIntervals.gr)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, BreakpointIntervals.gr, as.character=F)
  seq = seq[width(seq)>20]
  names(seq) = BreakpointIntervals.gr$ID
  writeXStringSet(seq, paste0("~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/", svcaller, "_OnlyNoPalmTreeTx_BreakpointSequences.fa"))

}


#### CIRCLE SEQ CIRCLE-CIRCLE CIRCLE-GENOME GENOME-GENOME CLASSIFICATION
WORKING_DIR = "~/Desktop/PalmTrees/Analysis/WorkspaceData/BreakpointSequences/"
generate_fa_file = function(this_tx, info, fname){
  BreakpointIntervals = data.frame(
    ID = c(paste0(this_tx$Sample, "_", this_tx$SVCaller, "_", info, "_",this_tx$ChrA, "_", this_tx$PosA), 
           paste0(this_tx$Sample, "_", this_tx$SVCaller, "_", info, "_",this_tx$ChrB, "_", this_tx$PosB)),
    chr = c(this_tx$ChrA, this_tx$ChrB),
    start = c(this_tx$PosA - ceiling(IntervalLength / 2), this_tx$PosB - ceiling(IntervalLength / 2)),
    end = c(this_tx$PosA + ceiling(IntervalLength / 2), this_tx$PosB + ceiling(IntervalLength / 2)))
  BreakpointIntervals$strand = "*"
  BreakpointIntervals = BreakpointIntervals %>% distinct()
  BreakpointIntervals.gr = makeGRangesFromDataFrame(df = BreakpointIntervals, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  BreakpointIntervals.gr = BreakpointIntervals.gr[!grepl("NoPalmTree", BreakpointIntervals.gr$ID)]
  seqlevels(BreakpointIntervals.gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
  seqinfo(BreakpointIntervals.gr) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
  BreakpointIntervals.gr <- trim(BreakpointIntervals.gr)
  seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, BreakpointIntervals.gr, as.character=F)
  seq = seq[width(seq)>20]
  names(seq) = BreakpointIntervals.gr$ID
  writeXStringSet(seq, fname)
}

IntervalLength = 60
svcallers = "Union" #unique(tx_circles$SVCaller)
circlemethods = c("CircleSeq", "WGSCircles") #unique(tx_circles$CircleMethod)
circlegenomeclasses = c("circle-circle", "circle-genome", "genome-genome")

for (svcaller in svcallers){
  for (circlemethod in circlemethods){
    for (circlegenomeclass in circlegenomeclasses){
      this_tx = tx_circles %>% 
        filter(CircleMethod == circlemethod, SVCaller == svcaller) %>%
        filter(grepl(",PRECISE", SVCallerInfo)) %>%
        filter(CircleGenomeClass == circlegenomeclass)
      generate_fa_file(this_tx, 
                       info = paste0(circlemethod, "_", circlegenomeclass),
                       fname = paste0(WORKING_DIR, svcaller, "_", circlemethod, "_", circlegenomeclass, "_BreakpointSequences.fa"))
      
      this_tx = tx_circles %>% 
        filter(CircleMethod == circlemethod, SVCaller == svcaller) %>%
        filter(grepl(",PRECISE", SVCallerInfo)) %>%
        filter(CircleGenomeClass != circlegenomeclass)
      generate_fa_file(this_tx, 
                       info = paste0(circlemethod, "_non-", circlegenomeclass),
                       fname = paste0(WORKING_DIR, svcaller, "_", circlemethod, "_non-", circlegenomeclass, "_BreakpointSequences.fa"))
    }
  }
}

