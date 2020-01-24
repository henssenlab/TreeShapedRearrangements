library(dplyr)

circles = read.table("/Volumes/Transcend/PalmTrees/Data/CircleseqCircles.bed")

circles = circles %>%
  dplyr::select(V1, V2, V3, V16, V17) 
colnames(circles) = c("chr", "start", "end", "sample", "class")
circles$length = abs(circles$end - circles$start)

circles %>%
  dplyr::select(-class) %>% 
  distinct() %>% 
  summarise(meanLength = mean(length, na.rm=T)) %>% 
  print()

circles %>%
  filter(class == "eccDNA") %>% 
  distinct() %>% 
  summarise(meanLength = mean(length, na.rm=T),
            sdLength  = sd(length, na.rm=T),
            medianLength = median(length, na.rm=T),
            n = n(),
            nSamples = n_distinct(sample)) %>% 
  print()

circles %>%
  filter(class == "ecDNA") %>% 
  distinct() %>% 
  summarise(meanLength = mean(length, na.rm=T),
            sdLength  = sd(length, na.rm=T),
            medianLength = median(length, na.rm=T),
            n = n(),
            nSamples = n_distinct(sample)) %>%
  print()

