library(dplyr)
library(ggplot2)
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan.Rdata")

pedpancan_meta = read.table("~/Desktop/PalmTrees/Data/PedPanCanMeta.csv", sep=";", header=T)
pedpancan_meta = pedpancan_meta %>% filter(seq_type == "wgs")

pedpancan_meta$hasPalmTree = pedpancan_meta$sample %in% unique(palmtrees$Sample)

mean(pedpancan_meta$hasPalmTree)

PalmTreesByEntity = pedpancan_meta %>%
  group_by(entity) %>%
  summarise(NumberOfSamples = n(), NumberOfPalmTreeSamples = sum(hasPalmTree), PercentPalmTreeSamples = 100 * NumberOfPalmTreeSamples / NumberOfSamples) %>%
  arrange(desc(PercentPalmTreeSamples)) %>%
  mutate(entity = factor(entity, levels=unique(entity)))

write.table(PalmTreesByEntity, file="~/Desktop/PalmTrees/Results/PedPanCan/PedPanCanPTNummbers.txt",
            quote = F,
            row.names = F, 
            sep = "\t")

PalmTreesByEntity %>%
  ggplot(aes(x=entity, y=PercentPalmTreeSamples, fill=entity)) + 
  geom_col() + 
  theme_kons1() +
  ylab("Samples with Palm Trees [%]") + 
  xlab("")+
  guides(fill=F) + 
  ggsave("~/Desktop/PalmTrees/Results/PedPanCan/PedPanCanPTNumbers.pdf", height = 3, width=4)

