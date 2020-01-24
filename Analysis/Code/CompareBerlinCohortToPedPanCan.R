rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree_PedPanCan.Rdata")

# How many neuroblastomas does the PedPanCan cohort have.
meta %>% filter(entity == "nb", seq_type=="wgs") %>% group_by(entity_sub) %>% summarise(n=n())

tx_original %>%
  group_by(Sample) %>% 
  summarise(ntx = n()) %>% 
  ungroup() %>% 
  left_join(meta %>% dplyr::rename(Sample = sample)) 

nrow(meta %>% filter(seq_type == "wgs"))
length(unique(tx_original$Sample))

npalmtrees = palmtrees %>%
  group_by(Sample) %>% 
  summarise(npt = n_distinct(PalmTreeID)) %>% 
  ungroup() 

all_samples = 
  meta %>% filter(seq_type == "wgs") %>% dplyr::rename(Sample = sample) %>% 
  dplyr::select(Sample, entity, entity_sub) %>% distinct() %>% 
  left_join(npalmtrees) %>%
  mutate(npt = ifelse(is.na(npt), 0, npt)) %>%
  mutate(entity = as.character(entity)) %>% 
  mutate(entity = ifelse(grepl("mb_", entity), "mb", entity)) %>%
  mutate(entity = ifelse(grepl("hgg", entity), "hgg", entity)) %>%
  mutate(entity = ifelse(grepl("b-all", entity), "b-all", entity)) %>%
  mutate(entity = ifelse(grepl("epd_", entity), "epd", entity)) %>%
  mutate(entity = toupper(entity))

all_samples %>%
  filter(entity=="NB") %>%
  group_by(entity_sub) %>%
  summarise(n = n(),
            nWithPT = sum(npt>0))

all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0)) %>% 
  mutate(entity = forcats::fct_reorder(entity, -percsamples_with_pt)) %>% 
  ggplot(aes(x=entity, y=percsamples_with_pt)) + 
  geom_col() + 
  theme_kons1()

all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0)) %>% 
  mutate(entity = forcats::fct_reorder(entity, -percsamples_with_pt)) %>% 
  ggplot(aes(x=entity, y=percsamples_with_pt)) + 
  geom_col() + 
  theme_kons1()

all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0)) %>% 
  gather("Class", "NSamples", 3:4) %>% 
  mutate(entity = forcats::fct_reorder(entity, -nsamples)) %>% 
  filter(nsamples > 0) %>% 
  mutate(entity = droplevels(entity)) %>% 
  mutate(Class = factor(Class, levels=c("nsamples_without_pt", "nsamples_with_pt"))) %>% 
  ggplot(aes(x=entity, y=NSamples, color=Class, fill=Class)) + 
  geom_col() + 
  theme_kons1() + 
  coord_flip()

all_samples %>%
  group_by(entity) %>%
  summarise(nsamples = n(),
            nsamples_with_pt = sum(npt > 0),
            nsamples_without_pt =  n()-sum(npt > 0),
            percsamples_with_pt = mean(npt > 0)) %>% 
  mutate(entity = forcats::fct_reorder(entity, percsamples_with_pt)) %>% 
  filter(nsamples > 0) %>% 
  mutate(entity = droplevels(entity)) %>% 
  ggplot(aes(x=entity, y=percsamples_with_pt)) + 
  geom_col() + 
  theme_kons1() + 
  coord_flip()

sum(all_samples$npt > 0) / 546
sum(all_samples$npt > 0) / 278

rm(list=ls())
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

n_samples_with_pt = length(unique(palmtrees$Sample))
n_samples_with_tx = length(unique(tx_original$Sample))
n_samples_with_pt/n_samples_with_tx

#### Compare our cohort to Pfister cohort

## 2 of 21 have PT in Pfister cohort
## 26 of 93 have PT in our cohort
prop.test(x=c(2,26), n=c(21,93))

#2 24
#21 72
fisher.test(matrix(data=c(2,24,21,72), nrow=2, ncol=2, byrow=T))

### Comparing our cohort (with clinical data) to the Pfister cohort ---- MYCN-amplfied

# 2 of 3 have PT in Pfister cohort
# 10 of 17 have PT in our cohort

prop.test(x=c(2,10), n=c(3,17))


#2 1
#10 7
fisher.test(matrix(data=c(2,1,10,7), nrow=2, ncol=2, byrow=T))
