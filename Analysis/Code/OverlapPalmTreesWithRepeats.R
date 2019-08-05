library(dplyr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg19)
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/Repeats.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")

tx_original = tx_original %>% as_tibble()
tx_original$RepeatA = NA
tx_original$RepeatB = NA
tx_original$RepeatClassA = NA
tx_original$RepeatClassB = NA
for (i in 1:nrow(tx_original)){
  print(i)
  repeat_a = repeats %>% filter(Chr == tx_original[[i,"ChrA"]],  tx_original[[i,"PosA"]]>=Start, tx_original[[i,"PosA"]]<=End) %>% .$Repeat %>% as.character() %>% as.list() %>% do.call(paste, .)
  repeat_b = repeats %>% filter(Chr == tx_original[[i,"ChrB"]],  tx_original[[i,"PosB"]]>=Start, tx_original[[i,"PosB"]]<=End) %>% .$Repeat %>% as.character() %>% as.list() %>% do.call(paste, .)
  repeat_class_a = repeats %>% filter(Chr == tx_original[[i,"ChrA"]],  tx_original[[i,"PosA"]]>=Start, tx_original[[i,"PosA"]]<=End) %>% .$Class %>% as.character() %>% as.list() %>% do.call(paste, .)
  repeat_class_b = repeats %>% filter(Chr == tx_original[[i,"ChrB"]],  tx_original[[i,"PosB"]]>=Start, tx_original[[i,"PosB"]]<=End) %>% .$Class %>% as.character() %>% as.list() %>% do.call(paste, .)
  if (length(repeat_a)>0){
    tx_original[[i, "RepeatA"]] = repeat_a
    tx_original[[i, "RepeatClassA"]] = repeat_class_a
  }
  if (length(repeat_b)>0){
    tx_original[[i, "RepeatB"]] = repeat_b
    tx_original[[i, "RepeatClassB"]] = repeat_class_b
  }}

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/tx_original_with_repeats.Rdata")

# # Generate random positions and get distribution
# random_positions = sample_random_regions(50000, 1) %>% mutate(Pos = Start) %>% select(Chr, Pos)
# random_positions$Repeat = NA
# random_positions$RepeatClass = NA
# for (i in 1:50){
#   print(i)
#   repeat_random_pos = repeats %>% filter(Chr == random_positions[i, "Chr"],  random_positions[i, "Pos"]>=Start, random_positions[i, "Pos"]<=End) %>% .$Repeat %>% as.character() %>% as.list() %>% do.call(paste, .)
#   repeat_class_random_pos = repeats %>% filter(Chr == random_positions[i, "Chr"],  random_positions[i, "Pos"]>=Start, random_positions[i, "Pos"]<=End) %>% .$Class %>% as.character() %>% as.list() %>% do.call(paste, .)
#   if (length(repeat_random_pos)>0){
#     random_positions[i, "Repeat"] = repeat_random_pos
#     random_positions[i, "RepeatClass"] = repeat_class_random_pos
#   }
# }

relative = repeats %>% mutate(Length = End-Start) %>% group_by(Repeat, Class) %>% summarise(OverallLength = sum(Length))
genome_length = sum(
  seqlengths(BSgenome.Hsapiens.UCSC.hg19)[standardChromosomes(BSgenome.Hsapiens.UCSC.hg19)]
) # this also counts Ns which is probably not correct
relative$ratio = relative$OverallLength / genome_length

relative$ratio_nexc = relative$OverallLength / 2897310462 # based on the number of non-N bases in hg19 according to http://genomewiki.ucsc.edu/index.php/Hg19_Genome_size_statistics
