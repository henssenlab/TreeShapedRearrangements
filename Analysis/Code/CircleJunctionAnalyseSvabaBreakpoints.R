rm(list=ls())
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Biostrings)
source("/Volumes/Transcend/PalmTrees/Analysis/Code/CustomThemes.R")
load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/SvabaCircleJunctions.Rdata")
load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/Circles.Rdata")

samples = unique(svaba_junctions$Sample)
for (sample in samples){
  svaba_junctions %>% 
    filter(Sample == sample) %>% 
    dplyr::select(ChrA,PosA,PosB) %>% 
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/SvabaCircles/", sample, "_SvabaPreciseJunctions.bed"),
                quote=F, sep="\t", col.names=F, row.names=F)
}

########################################################
##########        SvabaCircleStats          ############
########################################################

svaba_junctions %>% 
  group_by(Sample) %>% 
  summarise(nCircles = n()) %>% 
  ungroup() %>% 
  mutate(Sample = forcats::fct_reorder(Sample,-nCircles)) %>% 
  ggplot(aes(x=Sample, y=nCircles))+
  geom_col(fill="grey25") +
  xlab("") +
  ylab("Number of Circles") +
  ggtitle("Svaba Circles") +
  theme_kons1() +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/SvabaCircles_NumberBySample.pdf", width=4, height=2, useDingbats=F)

# Correlate Circle-seq vs Svaba Call call number in general
CircleSeqNumberOfCircles = circles %>%
  filter(Method=="CircleSeq") %>% 
  group_by(Sample) %>% 
  summarise(CircleSeqCircles = n()) %>% 
  ungroup()
svaba_junctions %>% 
  group_by(Sample) %>% 
  summarise(SvabaCircles = n()) %>% 
  ungroup() %>% 
  full_join(CircleSeqNumberOfCircles) %>% 
  ggplot(aes(x=CircleSeqCircles, y=SvabaCircles)) +
  geom_point() +
  geom_text_repel(aes(label=Sample), size=3) +
  theme_kons1() + 
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/SvabaCircles_NumberBySample_CorrelationWithCircleSeq.pdf", width=5, height=5, useDingbats=F)


svaba_junctions %>%
  mutate(CircleLength = PosB-PosA) %>% 
  ggplot(aes(x=CircleLength)) +
  geom_density(fill="grey25", color=NA) + 
  scale_x_continuous(trans="log1p", breaks=c(100,300,1000,3000,10000)) +
  theme_kons1() +
  xlab("Circle Length [bp]") +
  ylab("") +
  ggtitle("Svaba Circle Lengths") +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/SvabaCircles_CircleLength.pdf", width=3, height=2, useDingbats=F)

svaba_junctions %>%
  mutate(CircleLength = PosB-PosA) %>% 
  ggplot(aes(x=CircleLength)) +
  geom_density(fill="grey25", color=NA) + 
  scale_x_continuous(trans="log1p", breaks=c(100,300,1000,3000,10000)) +
  theme_kons1() +
  xlab("Circle Length [bp]") +
  ylab("") +
  ggtitle("Svaba Circle Lengths") +
  facet_wrap(Sample~.) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/SvabaCircles_CircleLength_bySample.pdf", width=10, height=10, useDingbats=F)


svaba_junctions %>%
  mutate(CircleLength = PosB-PosA) %>% 
  filter(CircleLength < 1500) %>% 
  ggplot(aes(x=CircleLength)) +
  geom_density(bw=3, fill="grey25", color=NA) + 
  theme_kons1() +
  xlab("Circle Length [bp]") +
  ylab("") +
  ggtitle("Svaba Circle Lengths") +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/SvabaCircles_CircleLength_Less1kb.pdf", width=3, height=2, useDingbats=F)

# n circles
svaba_junctions %>%
  mutate(CircleLength = PosB-PosA) %>% 
  filter(CircleLength < 1500) %>%
  nrow()

# how many samples
svaba_junctions$Sample %>% unique %>% length
  
##################################################
##########        HOMOLOGIES          ############
##################################################

svaba_junctions %>%
  mutate(hasHomology = HomologyLength>0) %>% 
  .$hasHomology %>% mean()
# 61.5% of precisely reconstructable circle junctions have an homology

svaba_junctions %>%
  mutate(hasHomology = HomologyLength>=5) %>% 
  .$hasHomology %>% mean()
# 6% of precisely reconstructable circle junctions have a homology of at least 5bp

svaba_junctions %>%
  mutate(hasHomology = HomologyLength>=10) %>% 
  .$hasHomology %>% mean()
# 0.4% of precisely reconstructable circle junctions have a homology of at least 5bp

svaba_junctions %>%
  ggplot(aes(x=HomologyLength)) +
  geom_histogram(fill="steelblue", color=NA) + 
  theme_kons1() +
  xlab("Homology Length at \n Circle Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/HomologyLength.pdf", width=3, height=2, useDingbats=F)

# n for paper
svaba_junctions %>% nrow()

svaba_junctions %>%
  ggplot(aes(x=HomologyLength)) +
  geom_histogram(fill="steelblue", color=NA) + 
  theme_kons1() +
  xlab("Homology Length at Circle Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  facet_wrap(Sample~.) + 
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/HomologyLength_bySample.pdf", width=6, height=6, useDingbats=F)

svaba_junctions %>%
  ggplot(aes(x=HomologyLength)) +
  geom_density(bw=2, fill="steelblue", color=NA) + 
  theme_kons1() +
  xlab("Homology Length at Circle Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  facet_wrap(Sample~.) + 
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/HomologyLength_bySample.pdf", width=6, height=6, useDingbats=F)

# Analyse Sequence
library(Biostrings)
svaba_junctions_seqs = svaba_junctions %>%
  dplyr::select(Sample, BNDPairID, HomologyLength, Homology) %>%
  distinct() %>%
  filter(HomologyLength>=5) %>% 
  dplyr::select(BNDPairID, Homology)
svaba_junctions_seqs_xstr = DNAStringSet(svaba_junctions_seqs$Homology)
names(svaba_junctions_seqs_xstr) = svaba_junctions_seqs$BNDPairID
writeXStringSet(svaba_junctions_seqs_xstr, paste0("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Homology_Min5.fa"))

## Analyse how that compares to randomized breakpoints
source("~/Desktop/PalmTrees/Analysis/Code/myHomology.R")
library(parallel)
cores=detectCores()

# This is computed on the cluster because I wrote insanely slow code ------------------------------------------------------------------

# svaba_junctions$myHomology = mcmapply(getHomologySvabaBreakpoints,
#                                        svaba_junctions$ChrA, svaba_junctions$PosA, svaba_junctions$DirectionA,
#                                        svaba_junctions$ChrB, svaba_junctions$PosB, svaba_junctions$DirectionB,
#                                        mc.cores=cores)
# svaba_junctions$myHomologyLength = nchar(svaba_junctions$myHomology)
# 
# unpermuted_svaba_junctions = svaba_junctions %>% filter(isPrecise) %>% dplyr::select(Sample, ChrA, PosA, DirectionA, ChrB, PosB, DirectionB, CircleGenomeClass, isInUnionPT) %>% distinct()
# 
# permuted_svaba_junctions = unpermuted_svaba_junctions
# B_Perm_i = sample.int(nrow(unpermuted_svaba_junctions),nrow(unpermuted_svaba_junctions),replace=F)
# permuted_svaba_junctions[, "ChrB"] = permuted_svaba_junctions[B_Perm_i, "ChrB"]
# permuted_svaba_junctions[, "PosB"] = permuted_svaba_junctions[B_Perm_i, "PosB"]
# permuted_svaba_junctions[, "DirectionB"] = permuted_svaba_junctions[B_Perm_i, "DirectionB"]
# permuted_svaba_junctions$myHomology = mcmapply(getHomologySvabaBreakpoints,
#                                                permuted_svaba_junctions$ChrA, permuted_svaba_junctions$PosA, permuted_svaba_junctions$DirectionA,
#                                                permuted_svaba_junctions$ChrB, permuted_svaba_junctions$PosB, permuted_svaba_junctions$DirectionB,
#                                                mc.cores=cores)
# permuted_svaba_junctions$myHomologyLength = nchar(permuted_svaba_junctions$myHomology)
# permuted_svaba_junctions$isPermuted = "Randomized"
# 
# unpermuted_svaba_junctions$myHomology = mcmapply(getHomologySvabaBreakpoints, 
#                                                  unpermuted_svaba_junctions$ChrA, unpermuted_svaba_junctions$PosA, unpermuted_svaba_junctions$DirectionA,
#                                                  unpermuted_svaba_junctions$ChrB, unpermuted_svaba_junctions$PosB, unpermuted_svaba_junctions$DirectionB,
#                                                  mc.cores=cores)
# unpermuted_svaba_junctions$myHomologyLength = nchar(unpermuted_svaba_junctions$myHomology)
# unpermuted_svaba_junctions$isPermuted = "Real"
# permuted_and_unpermuted_svaba_junctions = bind_rows(permuted_svaba_junctions, unpermuted_svaba_junctions)

# Load the results from the cluster
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/CircleJucntionAnalyseHomologyLengthVsRandom_OnCluster.Rdata")

ggplot(data=svaba_junctions, aes(x=HomologyLength, y=myHomologyLength)) +
  geom_point(size=1) +
  xlab(paste0("Svaba Homology Length \n\n",
              "Percent Sequence Identical: ", sprintf("%.2f",mean(svaba_junctions$Homology==svaba_junctions$myHomology, na.rm=T)), "\n",
              "Percent Length Identical: ", sprintf("%.2f",mean(svaba_junctions$HomologyLength==svaba_junctions$myHomologyLength, na.rm=T)), "\n",
              "SvabaHomology Percent >0bp: ", sprintf("%.2f",mean(svaba_junctions$HomologyLength>0, na.rm=T)), "\n",
              "MyHomology Percent >0bp: ", sprintf("%.2f",mean(svaba_junctions$myHomologyLength>0, na.rm=T)), "\n",
              "SvabaHomology Percent >=5bp: ", sprintf("%.2f",mean(svaba_junctions$HomologyLength>=5, na.rm=T)), "\n",
              "MyHomology Percent >=5bp: ", sprintf("%.2f",mean(svaba_junctions$myHomologyLength>=5, na.rm=T)), "\n",
              "SvabaHomology Percent >=10bp: ", sprintf("%.2f",mean(svaba_junctions$HomologyLength>=10, na.rm=T)), "\n",
              "MyHomology Percent >=10bp: ", sprintf("%.2f",mean(svaba_junctions$myHomologyLength>=10, na.rm=T)), "\n",
              "Mean SvabaHomology: ", sprintf("%.2f",mean(svaba_junctions$HomologyLength, na.rm=T)), "\n",
              "Mean myHomology: ", sprintf("%.2f",mean(svaba_junctions$myHomologyLength, na.rm=T)), "\n",
              "SD SvabaHomology: ",sprintf("%.2f",sd(svaba_junctions$HomologyLength, na.rm=T)), "\n",
              "SD myHomology: ", sprintf("%.2f",sd(svaba_junctions$myHomologyLength, na.rm=T)), "\n",
              "Median SvabaHomology: ", sprintf("%.2f",median(svaba_junctions$HomologyLength, na.rm=T)), "\n",
              "Median myHomology: ", sprintf("%.2f",median(svaba_junctions$myHomologyLength, na.rm=T)), "\n",
              "Mean Abs. Difference:", sprintf("%.2f", mean(abs(svaba_junctions$HomologyLength-svaba_junctions$myHomologyLength), na.rm=T)), "\n",
              "Median Abs. Difference:", sprintf("%.2f", median(abs(svaba_junctions$HomologyLength-svaba_junctions$myHomologyLength), na.rm=T)), "\n",
              "Correlation: ", sprintf("%.2f",cor.test(svaba_junctions$HomologyLength, svaba_junctions$myHomologyLength, use="complete.obs")$estimate)
  )) +
  ggtitle("Svaba Circle Junctions") +
  theme_kons1() +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/myHomologyVsSvabaHomologyScatter.pdf", height=8, width=5.5)

ggplot(data=svaba_junctions, aes(x=HomologyLength-myHomologyLength)) +
  geom_histogram(binwidth = 1) +
  ggtitle("Svaba Circle Junctions") +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/myHomologyVsSvabaHomology.pdf")

permuted_and_unpermuted_svaba_junctions %>%
  group_by(isPermuted) %>%
  summarise(MedianHomologyLength = median(myHomologyLength, na.rm=T),
            MeanHomologyLength = mean(myHomologyLength, na.rm=T))

t.test(myHomologyLength~isPermuted, data=permuted_and_unpermuted_svaba_junctions, alternative="two.sided") %>%
  capture.output() %>% 
  writeLines("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/HomologiesVsRandom_ttest.txt")

permuted_and_unpermuted_svaba_junctions %>% 
  ggplot(aes(x=myHomologyLength, fill=isPermuted)) + 
  geom_density(bw=2, alpha=0.5) +
  xlab("Homology Lengths") +
  ylab("Density") +
  xlim(0,30) +
  theme_kons1() +
  ggtitle("Svaba Circle Junctions") + 
  scale_color_manual(values=c("Randomized"="grey50", "Real"="red")) +
  scale_fill_manual(values=c("Randomized"="grey50", "Real"="red")) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/HomologiesVsRandomDensity.pdf", height=2.5, width=4)

permuted_and_unpermuted_svaba_junctions %>% 
  ggplot(aes(x=isPermuted, y=myHomologyLength, fill=isPermuted,color=isPermuted)) + 
  geom_jitter(alpha=0.1, size=1) +
  xlab("Homology Lengths") +
  ylab("Density") +
  theme_kons1() +
  scale_y_continuous(trans="log1p") +
  scale_color_manual(values=c("Randomized"="grey50", "Real"="red")) +
  scale_fill_manual(values=c("Randomized"="grey50", "Real"="red")) +
  ggtitle("Svaba Circle Junctions") + 
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/HomologiesVsRandomDots.pdf", height=4, width=5)


##################################################
##########        INSERTIONS          ############
##################################################

svaba_junctions %>%
  mutate(hasInsertion = InsertionLength>0) %>% 
  .$hasInsertion %>% mean()
# 2.8% of precisely reconstructable circle junctions have an insertion

svaba_junctions %>%
  mutate(hasInsertion = InsertionLength>=5) %>% 
  .$hasInsertion %>% mean()
# 1.8% of precisely reconstructable circle junctions have an insertions of at least 5bp

svaba_junctions %>%
  mutate(hasInsertion = InsertionLength>=10) %>% 
  .$hasInsertion %>% mean()
# 0.8% of precisely reconstructable circle junctions have an insertions of at least 5bp

svaba_junctions %>%
  distinct() %>%
  .$InsertionLength %>% table()

svaba_junctions %>%
  distinct() %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_histogram(fill="firebrick2", color=NA) + 
  theme_kons1() +
  xlab("Insertion Length at Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/InsertionLengths.pdf", height=2, width=4)

svaba_junctions %>%
  distinct() %>% 
  filter(InsertionLength>0) %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_histogram(fill="firebrick2", color=NA) + 
  theme_kons1() +
  xlab("Insertion Length at \n Circle Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/InsertionLengthsLargerZero.pdf", height=2, width=3)


svaba_junctions %>%
  distinct() %>% 
  filter(InsertionLength>0) %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_density(bw=1,fill="firebrick2", color=NA) + 
  theme_kons1() +
  xlab("Insertion Length at \n Circle Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/InsertionLengthsLargerZero_Density.pdf", height=2, width=3)

# n to report in the paper
svaba_junctions %>%
  distinct() %>% 
  filter(InsertionLength>0) %>%
  nrow()


svaba_junctions %>%
  distinct() %>% 
  filter(InsertionLength>0) %>% 
  ggplot(aes(x=InsertionLength, color=Sample)) +
  geom_density(bw=1,fill=NA) + 
  scale_color_grey() +
  theme_kons1() +
  xlab("Insertion Length at \n Circle Junction [bp]") +
  ylab("Count") +
  ggtitle("Svaba Circle Junctions") +
  guides(color=F, fill=F)+
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/InsertionLengthsLargerZero_bySample_Density.pdf", height=2, width=3)

library(Biostrings)
svaba_junctions_seqs = svaba_junctions %>%
  dplyr::select(Sample, ChrA, PosA, ChrB, PosB, InsertionLength, Insertion) %>%
  mutate(ID = paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB)) %>% 
  distinct() %>%
  filter(InsertionLength>=20) %>% 
  dplyr::select(ID, Insertion) 
svaba_junctions_seqs_xstr = DNAStringSet(svaba_junctions_seqs$Insertion)
names(svaba_junctions_seqs_xstr) =svaba_junctions_seqs$ID
writeXStringSet(svaba_junctions_seqs_xstr, paste0("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Insertions_Min20bp.fa"))

svaba_junctions_seqs = svaba_junctions %>%
  dplyr::select(Sample, ChrA, PosA, ChrB, PosB, InsertionLength, Insertion) %>%
  mutate(ID = paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB)) %>% 
  distinct() %>%
  filter(InsertionLength>=10) %>% 
  dplyr::select(ID, Insertion) 
svaba_junctions_seqs_xstr = DNAStringSet(svaba_junctions_seqs$Insertion)
names(svaba_junctions_seqs_xstr) =svaba_junctions_seqs$ID
writeXStringSet(svaba_junctions_seqs_xstr, paste0("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Insertions_Min10bp.fa"))

#system("source activate motifs; meme ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Insertions_Min10bp.fa -revcomp -nmotifs 10 -dna -oc ~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_Insertions_Min10bp_MEME")


#############################################################
##########        MOTIFS AT BREAKPOINTS          ############
#############################################################

IntervalLength = 20

svaba_junctions = svaba_junctions %>% 
  mutate(ID = paste0(Sample, "_", ChrA, ":", PosA, "_", ChrB, ":", PosB))

BreakpointIntervals = data.frame(
  ID = c(paste0(svaba_junctions$ID, "_Start"), paste0(svaba_junctions$ID, "_End")),
  chr = c(svaba_junctions$ChrA, svaba_junctions$ChrB),
  start = c(svaba_junctions$PosA - ceiling(IntervalLength / 2), svaba_junctions$PosB - ceiling(IntervalLength / 2)),
  end = c(svaba_junctions$PosA + ceiling(IntervalLength / 2), svaba_junctions$PosB + ceiling(IntervalLength / 2)))
BreakpointIntervals$strand = "*"
BreakpointIntervals = BreakpointIntervals %>% distinct()

BreakpointIntervals.gr = makeGRangesFromDataFrame(df = BreakpointIntervals, seqnames.field = "chr", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
seqlevels(BreakpointIntervals.gr, pruning.mode="coarse") = seqlevels(BSgenome.Hsapiens.UCSC.hg19)
seqinfo(BreakpointIntervals.gr) = seqinfo(BSgenome.Hsapiens.UCSC.hg19)
BreakpointIntervals.gr <- trim(BreakpointIntervals.gr)
seq = getSeq(BSgenome.Hsapiens.UCSC.hg19, BreakpointIntervals.gr, as.character=F)
seq = seq[width(seq)>20]
names(seq) = BreakpointIntervals.gr$ID
writeXStringSet(seq, paste0("~/Desktop/PalmTrees/Results/Figures/CircleJunctionAnalysis/CircleJunction_BreakpointIntervals_41bp.fa"))
