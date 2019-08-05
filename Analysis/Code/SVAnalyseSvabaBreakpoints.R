rm(list=ls())

library(dplyr)
library(tidyr)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")
source("~/Desktop/PalmTrees/Analysis/Code/ParseAllSV.R")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/SvabaSV.Rdata")
svaba_ins_hom = svaba_allsv %>% 
  select(Sample, BNDPairID, Insertion, InsertionLength, Homology, HomologyLength,
         ChrA, PosA, ChrB, PosB, DirectionA, DirectionB) %>%
  mutate(SampleBNDPairID = paste0(Sample, "_", BNDPairID))

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/MergedTx.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

union_txptinfo = txptinfo %>% filter(SVCaller == "Union")

svaba_wgscircles = 
  svaba_wgscircles %>% 
  mutate(isInUnionPT = BPID %in% union_txptinfo$BPID) %>% 
  select(Sample, BPID, isInUnionPT, CircleGenomeClass, CircleMethod, Precision, ID) %>%
  filter(Precision == "PRECISE") %>% 
  mutate(SampleBNDPairID = paste0(Sample, "_", ID)) %>% 
  left_join(svaba_ins_hom, by=c("Sample", "SampleBNDPairID"))

svaba_circleseq = 
  svaba_circleseq %>% 
  mutate(isInUnionPT = BPID %in% union_txptinfo$BPID) %>%
  select(Sample, BPID, isInUnionPT, CircleGenomeClass, CircleMethod, Precision, ID) %>%
  filter(Precision == "PRECISE") %>% 
  mutate(SampleBNDPairID = paste0(Sample, "_", ID)) %>% 
  left_join(svaba_ins_hom, by=c("Sample", "SampleBNDPairID"))

##################################################
##########        HOMOLOGIES          ############
##################################################


### General Microhomologies

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength) %>%
  distinct() %>%
  mutate(hasHomology = HomologyLength>0) %>% 
  .$hasHomology %>% mean()
# 68.0% havea microhomology

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength) %>%
  distinct() %>%
  mutate(hasHomology = HomologyLength>=5) %>% 
  .$hasHomology %>% mean()
# 10.0% have a microhomology of at least 5bp

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength) %>%
  distinct() %>%
  mutate(hasHomology = HomologyLength>=10) %>% 
  .$hasHomology %>% mean()
# 4.1% have a microhomology of at least 10bp

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength) %>%
  distinct() %>%
  .$HomologyLength %>% table()

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength) %>%
  distinct() %>% 
  ggplot(aes(x=HomologyLength)) +
  geom_density(fill="steelblue", color=NA) + 
  theme_kons1() + 
  xlab("Homology Length at Breakpoint [bp]") +
  ylab("Density")

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength) %>%
  distinct() %>% 
  ggplot(aes(x=HomologyLength)) +
  geom_histogram(fill="steelblue", color=NA) + 
  theme_kons1() +
  xlab("Homology Length at Breakpoint [bp]") +
  ylab("Count") +
  ggtitle("Svaba Interchr. Rearrangements") +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/HomologyLength.pdf", width=3, height=2, useDingbats=F)

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, HomologyLength) %>%
  distinct() %>% 
  ggplot(aes(x=HomologyLength, color=Sample)) +
  geom_density(bw=5,fill=NA) + 
  theme_kons1() +
  xlab("Homology Length at Breakpoint [bp]") +
  ylab("Count") +
  guides(color=F) +
  scale_color_grey() +
  ggtitle("Svaba Interchr. Rearrangements \n Homology Lengths By Sample") +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/HomologyLengthBySample.pdf", width=3, height=2, useDingbats=F)

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, HomologyLength) %>%
  distinct() %>% 
  group_by(Sample) %>% 
  summarise(meanHomologyLength = mean(HomologyLength, na.rm=T),
            nTx = n()) %>%
  arrange(desc(meanHomologyLength)) 

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, HomologyLength) %>%
  distinct() %>%
  arrange(desc(HomologyLength)) %>% 
  View()

### PT vs. Non Palm Tree
svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength, isInUnionPT) %>%
  distinct() %>% 
  ggplot(aes(y=HomologyLength, x=isInUnionPT, fill=isInUnionPT, color=isInUnionPT)) +
  geom_violin(bw=3, color=NA, alpha=0.5) + 
  geom_jitter(width=0.2) +
  theme_kons1() +
  ylab("Homology Length at Breakpoint [bp]") +
  scale_fill_manual(values=c("TRUE"="firebrick1", "FALSE"="steelblue")) +
  scale_color_manual(values=c("TRUE"="firebrick1", "FALSE"="steelblue")) +
  xlab("Breakpoint associated with a Cluster of Rearrangements") +
  guides(fill=F, color=F) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/UnionCalls_SvabaHomologies_UnionPTvsNonPT.pdf", height=4, width=4, useDingbats=F)

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength, isInUnionPT) %>%
  distinct() %>% 
  t.test(HomologyLength~isInUnionPT, data=., alternative="two.sided", paired=F) %>%
  print() %>% capture.output() %>% writeLines("~/Desktop/PalmTrees/Results/Figures/Homologies/UnionCalls_SvabaHomologies_UnionPTvsNonPT_ttest.txt")

### Homology Lengths -- WGSCircles CircleGenome Classification 
circle_genome_class_colors = c("genome-genome" = "steelblue", "circle-genome"="firebrick2", "circle-circle"="green4")

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength, CircleGenomeClass) %>%
  distinct() %>% 
  ggplot(aes(y=HomologyLength, x=CircleGenomeClass, fill=CircleGenomeClass, color=CircleGenomeClass)) +
  geom_violin(bw=3, alpha=0.5, draw_quantiles = 0.5) + 
  geom_jitter(width=0.2) +
  theme_kons1() +
  ylab("Homology Length at Breakpoint [bp]") +
  scale_fill_manual(values=circle_genome_class_colors) +
  scale_color_manual(values=circle_genome_class_colors) +
  xlab("") +
  guides(fill=F, color=F) +
  ggtitle("Microhomology Length (Svaba) by\nCircle/Genome Classification (WGS)")+
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/UnionCalls_SvabaHomologies_WGSCirclesCircleGenomeClass.pdf", height=4, width=4, useDingbats=F)

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength, CircleGenomeClass) %>%
  distinct() %>% 
  aov(HomologyLength~CircleGenomeClass, data=.) %>% 
  summary() %>% capture.output() %>% writeLines("~/Desktop/PalmTrees/Results/Figures/Homologies/UnionCalls_SvabaHomologies_WGSCirclesCircleGenomeClass_Anova.txt")

svaba_wgscircles %>%
  dplyr::select(BPID, HomologyLength, CircleGenomeClass) %>%
  distinct() %>% 
  aov(HomologyLength~CircleGenomeClass, data=.) %>% 
  TukeyHSD %>% capture.output() %>% writeLines("~/Desktop/PalmTrees/Results/Figures/Homologies/UnionCalls_SvabaHomologies_WGSCirclesCircleGenomeClass_Anova_TukeyHSD.txt")

### CircleSeq CircleGenome Classification

svaba_circleseq %>%
  dplyr::select(BPID, HomologyLength, CircleGenomeClass) %>%
  distinct() %>% 
  ggplot(aes(y=HomologyLength, x=CircleGenomeClass, fill=CircleGenomeClass, color=CircleGenomeClass)) +
  geom_violin(bw=3, alpha=0.5, draw_quantiles = 0.5) + 
  #stat_summary(fun.y=mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar") +
  geom_jitter(width=0.2) +
  theme_kons1() +
  ylab("Homology Length at Breakpoint [bp]") +
  scale_fill_manual(values=circle_genome_class_colors) +
  scale_color_manual(values=circle_genome_class_colors) +
  xlab("") +
  guides(fill=F, color=F) +
  ggtitle("Microhomology Length (Svaba) by\nCircle/Genome Classification (Circle-seq)")+
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/UnionCalls_SvabaHomologies_CircleSeqCirclesCircleGenomeClass.pdf", height=4, width=4, useDingbats=F)

# Analyse Sequence
library(Biostrings)
svaba_wgscircles_seqs = svaba_wgscircles %>%
  dplyr::select(Sample, BPID, HomologyLength, Homology) %>%
  distinct() %>%
  filter(HomologyLength>=5) %>% 
  dplyr::select(BPID, Homology) 
svaba_wgscircles_seqs_xstr = DNAStringSet(svaba_wgscircles_seqs$Homology)
names(svaba_wgscircles_seqs_xstr) =svaba_wgscircles_seqs$BPID
writeXStringSet(svaba_wgscircles_seqs_xstr, paste0("~/Desktop/PalmTrees/Results/Figures/Homologies/SvabaTx_Homologies.fa"))

## Analyse how that compares to randomized breakpoints
source("~/Desktop/PalmTrees/Analysis/Code/myHomology.R")
library(parallel)
cores=detectCores()

svaba_wgscircles$myHomology = mcmapply(getHomologySvabaBreakpoints,
                            svaba_wgscircles$ChrA, svaba_wgscircles$PosA, svaba_wgscircles$DirectionA,
                            svaba_wgscircles$ChrB, svaba_wgscircles$PosB, svaba_wgscircles$DirectionB,
                            mc.cores=cores)
svaba_wgscircles$myHomologyLength = nchar(svaba_wgscircles$myHomology)
ggplot(data=svaba_wgscircles, aes(x=HomologyLength, y=myHomologyLength)) +
  geom_point(size=1) +
  xlab(paste0("Svaba Homology Length \n\n",
             "Percent Sequence Identical: ", sprintf("%.2f",mean(svaba_wgscircles$Homology==svaba_wgscircles$myHomology, na.rm=T)), "\n",
             "Percent Length Identical: ", sprintf("%.2f",mean(svaba_wgscircles$HomologyLength==svaba_wgscircles$myHomologyLength, na.rm=T)), "\n",
             "SvabaHomology Percent >0bp: ", sprintf("%.2f",mean(svaba_wgscircles$HomologyLength>0, na.rm=T)), "\n",
             "MyHomology Percent >0bp: ", sprintf("%.2f",mean(svaba_wgscircles$myHomologyLength>0, na.rm=T)), "\n",
             "SvabaHomology Percent >=5bp: ", sprintf("%.2f",mean(svaba_wgscircles$HomologyLength>=5, na.rm=T)), "\n",
             "MyHomology Percent >=5bp: ", sprintf("%.2f",mean(svaba_wgscircles$myHomologyLength>=5, na.rm=T)), "\n",
             "SvabaHomology Percent >=10bp: ", sprintf("%.2f",mean(svaba_wgscircles$HomologyLength>=10, na.rm=T)), "\n",
             "MyHomology Percent >=10bp: ", sprintf("%.2f",mean(svaba_wgscircles$myHomologyLength>=10, na.rm=T)), "\n",
             "Mean SvabaHomology: ", sprintf("%.2f",mean(svaba_wgscircles$HomologyLength, na.rm=T)), "\n",
             "Mean myHomology: ", sprintf("%.2f",mean(svaba_wgscircles$myHomologyLength, na.rm=T)), "\n",
             "SD SvabaHomology: ",sprintf("%.2f",sd(svaba_wgscircles$HomologyLength, na.rm=T)), "\n",
             "SD myHomology: ", sprintf("%.2f",sd(svaba_wgscircles$myHomologyLength, na.rm=T)), "\n",
             "Median SvabaHomology: ", sprintf("%.2f",median(svaba_wgscircles$HomologyLength, na.rm=T)), "\n",
             "Median myHomology: ", sprintf("%.2f",median(svaba_wgscircles$myHomologyLength, na.rm=T)), "\n",
             "Mean Abs. Difference:", sprintf("%.2f", mean(abs(svaba_wgscircles$HomologyLength-svaba_wgscircles$myHomologyLength), na.rm=T)), "\n",
             "Median Abs. Difference:", sprintf("%.2f", median(abs(svaba_wgscircles$HomologyLength-svaba_wgscircles$myHomologyLength), na.rm=T)), "\n",
             "Correlation: ", sprintf("%.2f",cor.test(svaba_wgscircles$HomologyLength, svaba_wgscircles$myHomologyLength, use="complete.obs")$estimate)
           )) +
  theme_kons1() +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/myHomologyVsSvabaHomologyScatter.pdf", height=8, width=5.5)

ggplot(data=svaba_wgscircles, aes(x=HomologyLength-myHomologyLength)) +
  geom_histogram(binwidth = 1) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/myHomologyVsSvabaHomology.pdf")

unpermuted_svaba = svaba_wgscircles %>% filter(Precision == "PRECISE") %>% dplyr::select(Sample, ChrA, PosA, DirectionA, ChrB, PosB, DirectionB, CircleGenomeClass, isInUnionPT) %>% distinct()

permuted_svaba = unpermuted_svaba
B_Perm_i = sample.int(nrow(unpermuted_svaba),nrow(unpermuted_svaba),replace=F)
permuted_svaba[, "ChrB"] = permuted_svaba[B_Perm_i, "ChrB"]
permuted_svaba[, "PosB"] = permuted_svaba[B_Perm_i, "PosB"]
permuted_svaba[, "DirectionB"] = permuted_svaba[B_Perm_i, "DirectionB"]
permuted_svaba$myHomology = mcmapply(getHomologySvabaBreakpoints,
                                     permuted_svaba$ChrA, permuted_svaba$PosA, permuted_svaba$DirectionA,
                                     permuted_svaba$ChrB, permuted_svaba$PosB, permuted_svaba$DirectionB,
                                     mc.cores=cores)
permuted_svaba$myHomologyLength = nchar(permuted_svaba$myHomology)
permuted_svaba$isPermuted = "Randomized"

unpermuted_svaba$myHomology = mcmapply(getHomologySvabaBreakpoints, 
                                       unpermuted_svaba$ChrA, unpermuted_svaba$PosA, unpermuted_svaba$DirectionA,
                                       unpermuted_svaba$ChrB, unpermuted_svaba$PosB, unpermuted_svaba$DirectionB,
                                       mc.cores=cores)
unpermuted_svaba$myHomologyLength = nchar(unpermuted_svaba$myHomology)
unpermuted_svaba$isPermuted = "Real"

permuted_and_unpermuted_svaba = bind_rows(permuted_svaba, unpermuted_svaba)

permuted_and_unpermuted_svaba %>%
  group_by(isPermuted) %>% 
  summarise(MedianHomologyLength = median(myHomologyLength),
            MeanHomologyLength = mean(myHomologyLength))

t.test(myHomologyLength~isPermuted, data=permuted_and_unpermuted_svaba, alternative="two.sided") %>%
  capture.output() %>% 
  writeLines("~/Desktop/PalmTrees/Results/Figures/Homologies/HomologiesVsRandom_ttest.txt")

permuted_and_unpermuted_svaba %>% 
  ggplot(aes(x=myHomologyLength, fill=isPermuted)) + 
  geom_density(bw=2, alpha=0.5) +
  xlab("Homology Lengths") +
  ylab("Density") +
  xlim(0,50) +
  theme_kons1() +
  scale_color_manual(values=c("Randomized"="grey50", "Real"="red")) +
  scale_fill_manual(values=c("Randomized"="grey50", "Real"="red")) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/HomologiesVsRandomDensity.pdf", height=4, width=5)

permuted_and_unpermuted_svaba %>% 
  ggplot(aes(x=isPermuted, y=myHomologyLength, fill=isPermuted,color=isPermuted)) + 
  geom_jitter(alpha=0.1, size=1) +
  xlab("Homology Lengths") +
  ylab("Density") +
  theme_kons1() +
  scale_y_continuous(trans="log1p") +
  scale_color_manual(values=c("Randomized"="grey50", "Real"="red")) +
  scale_fill_manual(values=c("Randomized"="grey50", "Real"="red")) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Homologies/HomologiesVsRandomDots.pdf", height=4, width=5)

##################################################
##########        INSERTIONS          ############
##################################################

### General Insertion

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  mutate(hasInsertion = InsertionLength>0) %>% 
  .$hasInsertion %>% mean()
# 23% of precisely reconstructable Brekapoints have an insertions

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  mutate(hasInsertion = InsertionLength>=5) %>% 
  .$hasInsertion %>% mean()
# 14.5% of precisely reconstructable Brekapoints have an insertions of at least 5bp

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  mutate(hasInsertion = InsertionLength>=10) %>% 
  .$hasInsertion %>% mean()
# 9.8% of precisely reconstructable Brekapoints have an insertions of at least 5bp

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>%
  .$InsertionLength %>% table()

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_density(fill="firebrick2", color=NA) + 
  theme_kons1() + 
  xlab("Insertion Length at Breakpoint [bp]") +
  ylab("Density")

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_histogram(fill="firebrick2", color=NA) + 
  theme_kons1() +
  xlab("Insertion Length at Breakpoint [bp]") +
  ylab("Count") +
  ggtitle("Svaba Interchr. Rearrangements") +
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Insertions/SvabaBreakpoints_InsertionLengths.pdf", height=3, width=2)

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  filter(InsertionLength>0) %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_histogram(fill="firebrick2", color=NA) + 
  theme_kons1() +
  xlab("Insertion Length \n at SV [bp]") +
  ylab("Count") +
  ggtitle("Svaba Interchr. SV") +
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Insertions/SvabaBreakpoints_InsertionLengths_LargerZero.pdf", height=3, width=2)

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength) %>%
  distinct() %>% 
  ggplot(aes(x=InsertionLength)) +
  geom_density(fill="firebrick2", color=NA) + 
  theme_kons1() +
  xlab("Insertion Length \n at SV [bp]") +
  ylab("") +
  ggtitle("Svaba Interchr. SV") +
  scale_x_continuous(trans="log1p", breaks = c(0,3,10,30,100,300)) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/Insertions/SvabaBreakpoints_InsertionLengths_Density.pdf", height=3, width=2)

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, InsertionLength) %>%
  distinct() %>% 
  ggplot(aes(x=InsertionLength, color=Sample)) +
  geom_histogram(bw=3,fill=NA) + 
  theme_kons1() +
  xlab("Insertion Length \n at SV [bp]") +
  ylab("Count") +
  guides(color=F) +
  scale_color_grey() 

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, InsertionLength) %>%
  distinct() %>% 
  group_by(Sample) %>% 
  summarise(meanInsertionLength = mean(InsertionLength, na.rm=T),
            nTx = n()) %>%
  arrange(desc(meanInsertionLength))

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, InsertionLength, Insertion) %>%
  distinct() %>%
  arrange(desc(InsertionLength)) %>% 
  View()

svaba_wgscircles %>%
  dplyr::select(Sample, BPID, InsertionLength, Insertion) %>%
  distinct() %>%
  arrange(desc(InsertionLength)) %>% 
  write.table("~/Desktop/PalmTrees/Results/Figures/Insertions/SvabaTx_Insertions.txt",
              col.names=T, row.names=F, quote=F, sep="\t")

library(Biostrings)
svaba_wgscircles_seqs = svaba_wgscircles %>%
  dplyr::select(Sample, BPID, InsertionLength, Insertion) %>%
  distinct() %>%
  filter(InsertionLength>0) %>% 
  dplyr::select(BPID, Insertion) 
svaba_wgscircles_seqs_xstr = DNAStringSet(svaba_wgscircles_seqs$Insertion)
names(svaba_wgscircles_seqs_xstr) =svaba_wgscircles_seqs$BPID
writeXStringSet(svaba_wgscircles_seqs_xstr, paste0("~/Desktop/PalmTrees/Results/Figures/Insertions/SvabaTx_Insertions.fa"))

svaba_wgscircles_seqs = svaba_wgscircles %>%
  dplyr::select(Sample, BPID, InsertionLength, Insertion) %>%
  distinct() %>%
  filter(InsertionLength>=20) %>% 
  dplyr::select(BPID, Insertion) 
svaba_wgscircles_seqs_xstr = DNAStringSet(svaba_wgscircles_seqs$Insertion)
names(svaba_wgscircles_seqs_xstr) =svaba_wgscircles_seqs$BPID
writeXStringSet(svaba_wgscircles_seqs_xstr, paste0("~/Desktop/PalmTrees/Results/Figures/Insertions/SvabaTx_Insertions_Min20bp.fa"))

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength, CircleGenomeClass) %>%
  distinct() %>% 
  aov(InsertionLength~CircleGenomeClass, data=.) %>% 
  summary() %>% capture.output() %>% writeLines("~/Desktop/PalmTrees/Results/Figures/Insertions/UnionCalls_SvabaInsertions_WGSCirclesCircleGenomeClass_Anova.txt")

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength, CircleGenomeClass) %>%
  distinct() %>% 
  aov(InsertionLength~CircleGenomeClass, data=.) %>% 
  TukeyHSD %>% capture.output() %>% writeLines("~/Desktop/PalmTrees/Results/Figures/Insertions/UnionCalls_SvabaInsertions_WGSCirclesCircleGenomeClass_Anova_TukeyHSD.txt")

svaba_wgscircles %>%
  dplyr::select(BPID, InsertionLength, isInUnionPT) %>%
  distinct() %>% 
  t.test(InsertionLength~isInUnionPT, data=., alternative="two.sided", paired=F) %>%
  print() %>% capture.output() %>% writeLines("~/Desktop/PalmTrees/Results/Figures/Insertions/UnionCalls_SvabaInsertions_UnionPTvsNonPT_ttest.txt")

# TODO: Blast that -- anything interesting?
# FIMO: too short for that.

