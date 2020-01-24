rm(list=ls())
library(dplyr)
library(ggplot2)
load("/Volumes/Transcend/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

# Add Palm Tree Size (#bp between first and last rearrangement breakpoint)
palmtrees = palmtrees %>% mutate(Length = LastElement-FirstElement)

# some *few* translocations are part of two palm trees, 
# that can happen if two clusters are far apart from each 
# other and one rearrangement is just in the middle of the
# other cluster. We decided not to merge these. Now, make sure that 
# they do not affect overall statistics
tx_original = distinct(tx_original)
txptinfo = txptinfo %>% dplyr::select(-PalmTreeID) %>% distinct() 

# Beautify...
SVCallerColors = RColorBrewer::brewer.pal(7, "Set1")
names(SVCallerColors) = unique(tx_original$SVCaller)

# Ext Fig. 7a
tx_original %>% 
  group_by(Sample, SVCaller) %>%
  summarise(nTx = n_distinct(BPID)) %>%
  ungroup() %>% 
  ggplot(aes(x=SVCaller, y=nTx, fill=SVCaller)) +
  geom_violin() + 
  scale_y_log10() +
  theme_classic() + 
  ylab("Interchr. Rearrangements\nper Patient") + 
  xlab("")+
  scale_fill_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 14), 
        axis.title =  element_text(size = 14), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/NumberTxDistr.pdf", height=4, width=4)

# 7a counts
tx_original %>% 
  group_by(SVCaller) %>%
  summarise(n=n_distinct(Sample))

# 7a updated
tx_original %>% 
  group_by(Sample, SVCaller) %>%
  summarise(nTx = n_distinct(BPID)) %>%
  ungroup() %>% 
  ggplot(aes(x=SVCaller, y=nTx, fill=SVCaller)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_y_log10() +
  theme_classic() + 
  ylab("Interchr. Rearrangements\nper Patient") + 
  xlab("")+
  scale_fill_manual(values=SVCallerColors) +
  scale_color_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 14), 
        axis.title =  element_text(size = 14), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("/Volumes/Transcend/PalmTrees/Results/Figures/NumberTxDistrUpdate.pdf", height=4, width=4)

# Ext Figure 7b
palmtrees %>% 
  #filter(SVCaller %notin% c("Brass", "Novobreak")) %>% 
  group_by(Sample, SVCaller) %>%
  summarise(nTx = n_distinct(PalmTreeID)) %>%
  ungroup() %>% 
  ggplot(aes(x=SVCaller, y=nTx, fill=SVCaller)) +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) + 
  geom_jitter(color="black", size=1) + 
  scale_y_log10() +
  theme_classic() + 
  ylab("Tree-shape Rearrangments\nper Tumor") + 
  xlab("")+
  scale_fill_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 14), 
        axis.title =  element_text(size = 14), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("/Volumes/Transcend/PalmTrees/Results/Figures/NumberPalmTreesDistrUpdate.pdf", height=4, width=4,
         useDingbats=F)

# counts
palmtrees %>% 
  #filter(SVCaller %notin% c("Brass", "Novobreak")) %>% 
  group_by(Sample, SVCaller) %>%
  summarise(nTx = n_distinct(PalmTreeID)) %>%
  ungroup() %>% 
  group_by(SVCaller) %>%
  summarise(n=n())

# Ext Fig 7c  
palmtrees %>% 
  mutate(SVCaller = factor(palmtrees$SVCaller, levels=c("Novobreak", "Svaba", "Delly", "Smufin", "Brass", "AtLeastTwo", "Union"))) %>%
  group_by(SVCaller) %>%
  distinct() %>%
  summarise(nPT = n_distinct(PalmTreeID)) %>%
  ggplot(aes(y=nPT, x=SVCaller, fill=SVCaller)) + 
  geom_col() + 
  theme_classic() + 
  ylab("Number of Tree-shaped\nRearrangements") + 
  xlab("")+
  scale_fill_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 14), 
        axis.title =  element_text(size = 14), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/NumberPTbyCallers.pdf", height=4, width=4)


palmtrees %>% 
  mutate(SVCaller = factor(palmtrees$SVCaller, levels=c("Novobreak", "Svaba", "Delly", "Smufin", "Brass", "AtLeastTwo", "Union"))) %>%
  group_by(SVCaller) %>%
  distinct() %>%
  ggplot(aes(y=n, x=SVCaller, fill=SVCaller)) + 
  geom_violin() + 
  theme_classic() + 
  ylab("Translocations per Palm Tree") + 
  scale_fill_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 10), 
        axis.title =  element_text(face="bold", size = 10), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/TxPerPTbyCallers.pdf", height=4, width=4)

palmtrees %>% 
  mutate(SVCaller = factor(palmtrees$SVCaller, levels=c("Novobreak", "Svaba", "Delly", "Smufin", "Brass", "AtLeastTwo", "Union"))) %>%
  group_by(SVCaller) %>%
  distinct() %>%
  ggplot(aes(y=LastElement-FirstElement, x=SVCaller, fill=SVCaller)) + 
  geom_violin() + 
  theme_classic() + 
  ylab("Length") + 
  scale_y_log10(labels=scales::comma) +
  scale_fill_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 10), 
        axis.title =  element_text(face="bold", size = 10), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/PTLengthbyCallers.pdf", height=4, width=4)

tx_original %>%
  group_by(SVCaller) %>%
  summarise(ntx = n_distinct(BPID)) %>%
  ggplot(aes(x=SVCaller, y=ntx, fill=SVCaller)) +
  geom_col() +
  theme_classic() + 
  ylab("# Translocations") + 
  scale_fill_manual(values=SVCallerColors) +
  theme(text = element_text(family = "Helvetica"),
        axis.text.y = element_text(color="black"),
        axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
        axis.text = element_text(size = 10), 
        axis.title =  element_text(face="bold", size = 10), 
        plot.title = element_text(face="plain", size = 12, hjust=0.5))+
  guides(fill=F) +
  ggsave("~/Desktop/PalmTrees/Results/Figures/NumberTxbyCallers.pdf", height=4, width=4)

ggplot(data=palmtrees, aes(x=Length)) + 
  geom_histogram(fill = "steelblue", color=NA) + 
  scale_x_log10(labels = scales::comma) +
  theme_classic() + 
  guides(color=FALSE, fill=FALSE) + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") + 
  xlab("Palm Tree Region Size") + 
  ylab("Number of Palm Trees")
ggsave("~/Desktop/PalmTrees/Results/Figures/PalmTreeLengthDistribution.pdf")


## that still includes only defining (after filtering before palm tree calling) but not all translocations
ggplot(data=palmtrees, aes(x=n)) + 
  geom_histogram(fill = "steelblue", color=NA) + 
  scale_x_continuous(labels = scales::comma) +
  theme_classic() + 
  guides(color=FALSE, fill=FALSE) + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") + 
  xlab("Number of translocations per Palm Tree") + 
  ylab("Count")
ggsave("~/Desktop/PalmTrees/Results/Figures/NumberTxPerPalmTree.pdf")

ggplot(data=palmtrees %>% group_by(Sample) %>% summarise(PTperSample = n_distinct(PalmTreeID)), aes(x=PTperSample)) + 
  geom_histogram(fill = "steelblue", color=NA) + 
  scale_x_continuous(labels = scales::comma) +
  theme_classic() + 
  guides(color=FALSE, fill=FALSE) + 
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") + 
  xlab("Number of Palm Trees per Sample") + 
  ylab("Count")
ggsave("~/Desktop/PalmTrees/Results/Figures/NumberPTperPatient.pdf")
