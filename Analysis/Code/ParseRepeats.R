rm(list=ls())
library(dplyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)

fnames = paste0("~/Desktop/PalmTrees/Data/UCSCRepetitveRegions/chromOut/",
       gsub("chr", "", standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)), "/",
       standardChromosomes(BSgenome.Hsapiens.UCSC.hg19), ".fa.out")

repeats = list()
for (i in 1:length(fnames)){
  repeats[[i]] = read.table(fnames[i], skip=3) %>% as_tibble()
}
repeats = do.call(rbind, repeats)
repeats = repeats %>% 
  mutate(Chr = V5,
         Start = V6,
         End = V7,
         Strand = V9,
         Repeat = V10,
         Class = V11) %>%
  select(Chr, Start, End, Repeat, Class)

repeats_gr = makeGRangesFromDataFrame(repeats, seqnames.field="Chr", start.field="Start", end.field="End", keep.extra.columns = T)

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/Repeats.Rdata")