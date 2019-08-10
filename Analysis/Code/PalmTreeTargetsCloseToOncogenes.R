library(dplyr)
library(ggplot2)
library(parallel)

# Load Palm Tree Data
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")

# Load COSMIC genes 
cosmic = read.table("~/Desktop/PalmTrees/Data/cancer_gene_census.csv", header=T, sep=",")

rm(list=setdiff(ls(), c("txptinfo", "cosmic")))
source("~/Desktop/PalmTrees/Analysis/Code/ParseSVCallerData.R")
source("~/Desktop/PalmTrees/Analysis/Code/RegionSampling.R")
source("~/Desktop/PalmTrees/Analysis/Code/CustomThemes.R")
cores <- detectCores()

# Prepare COSMIC table
cosmic$Chr = paste0("chr", gsub(":.*", "", cosmic$Genome.Location))
cosmic$Start = cosmic$Genome.Location %>% gsub("^.*:", "", .) %>% gsub("-.*", "", .) %>% as.numeric()
cosmic$End = gsub("^.*-", "", cosmic$Genome.Location) %>% as.numeric()


# Analysis is performed separately on three gene sets (geneclasses)
# - all COSMIC genes
# - COSMIC Oncogenes
# - COSMIC Tumor Suppressor Genes
cosmic_oncogenes = cosmic %>% filter(!is.na(Start), !is.na(End),
                                     isdefchrom2(Chr),
                                     grepl("oncogene", cosmic$Role.in.Cancer))
cosmic_tsgs = cosmic %>% filter(!is.na(Start), !is.na(End),
                                isdefchrom2(Chr),
                                grepl("TSG", cosmic$Role.in.Cancer))
cosmic_all = cosmic %>% filter(!is.na(Start), !is.na(End),
                               isdefchrom2(Chr))
geneclasses = list(cosmic_oncogenes = cosmic_oncogenes,
                   cosmic_tsgs = cosmic_tsgs,
                   cosmic_all = cosmic_all)


# Function dist_to_closest_gene returns the distance of Position 
# "Chr":"Pos" to closest gene in table "this_cosmic" and returns
# -1 whenever a breakpoint disrupts (= locates within) a gene body
dist_to_closest_gene = function(Chr, Pos, this_cosmic){
  
  # Make sure input behaves well
  if (length(Chr)!=1) stop("length(Chr) != 1")
  if (length(Pos)!=1) stop("length(Pos) != 1")
  if (sum(this_cosmic$End < this_cosmic$Start) > 0) stop("Cosmic Coordinates: At least one End < Start")
  if (!isdefchrom(Chr)) stop("Non-standard input chromosome")
  
  thischr = Chr
  cosmic = this_cosmic %>%
    filter(Chr == thischr) %>%
    mutate(DistStart = (Start-Pos),
           DistEnd = (Pos-End),
           isDisrupted = ((DistStart < 0) & (DistEnd < 0)))
  
  if (sum(cosmic$isDisrupted, na.rm=T) > 0) return(-1)
  diststarts = cosmic$DistStart
  distends = cosmic$DistEnd
  smallest_positive_dist = min(c(diststarts[diststarts>0], distends[distends>0]))
  
  return(smallest_positive_dist)
}

# #Testing...
# dist_to_closest_gene("chr2", 16084000, cosmic_oncogenes) == -1 # disrupts MYCN, should return -1
# dist_to_closest_gene("chr2", 16082182, cosmic_oncogenes) == 5 # distance to MYCN (coding region) start is 5
# dist_to_closest_gene("chr2", 16086319, cosmic_oncogenes) == 100 # distance to MYCN (coding region) end is 100


# ------------------------------------------------------------------
# -------------- CREATE RANDOMIZED BREAKPOINTS ---------------------
# ------------------------------------------------------------------


# Pre-allocating
percent_disrupt = list()
distance_distribution = list()
mean_distance = list()
sd_distance = list()
median_distance = list()
random_positions = list()
random_positions_percent_disrupt = list()
random_positions_median_distance = list()


# Perform Randomization separately for all COSCIC genes vs. 
# Oncogenes only vs. TSGs only
for (gci in 1:length(geneclasses)){
  
  
  # For actual palm tree associated brekapoints, calculate distance to closest 
  # COSMIC gene and whether brekapoints disrupts a COSMIC gene
  txptinfo$DistToClosestOncogene = mapply(function (c, p) dist_to_closest_gene(c,p,geneclasses[[gci]]), txptinfo$TargetChrom, txptinfo$TargetPos)
  txptinfo$DisruptsCosmicOncogene = (txptinfo$DistToClosestOncogene == -1)
  
  
  # Consider the merged Tx Dataset "Union"only
  svcallers = c("Union")
  
  
  # Pre-allocating
  percent_disrupt[[gci]] = list()
  distance_distribution[[gci]] = list()
  mean_distance[[gci]] = list()
  sd_distance[[gci]] = list()
  median_distance[[gci]] = list()
  random_positions[[gci]] = list()
  random_positions_percent_disrupt[[gci]] = list()
  random_positions_median_distance[[gci]] = list()
  
  
  # If we would anaylse different Tx call sets, we could loop through them
  for (svi in 1:length(svcallers)){
    
    
    # Calculate summary statistics for real data
    this_txpt = txptinfo %>% filter(SVCaller == svcallers[svi])
    percent_disrupt[[gci]][[svi]] = mean(this_txpt$DistToClosestOncogene == -1)
    distance_distribution[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1)
    mean_distance[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>% .$DistToClosestOncogene %>% mean()
    sd_distance[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>% .$DistToClosestOncogene %>% sd()
    median_distance[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>% .$DistToClosestOncogene %>% median()
    
    
    # Sample 500 random, synthetic datasets from whitelisted genome and compute
    # for all breakpoints the distance to the closest COSMIC gene and wheter
    # it disrupts a COSMIC gene
    reps = 500
    this_random_positions = sample_random_regions_broad(n=reps*nrow(this_txpt), segment_length=0)
    this_random_positions$DistToClosestOncogene = mclapply(1:nrow(this_random_positions), function (i) dist_to_closest_gene(this_random_positions[[i,"Chr"]],this_random_positions[[i, "Start"]],geneclasses[[gci]]), mc.cores=cores) %>% as.numeric()
    this_random_positions$Rep = rep(1:reps, nrow(this_txpt)) %>% as.factor() # ID of Randomized Dataset
    this_random_positions$isDisruption = (this_random_positions$DistToClosestOncogene == -1)
    this_random_positions$isRandomDraw = TRUE
    
    
    # Store all randomized data
    random_positions[[gci]][[svi]] = this_random_positions
    
    
    # Store summary statistics for each of the 500 randomized datasets
    random_positions_percent_disrupt[[gci]][[svi]] = this_random_positions %>% group_by(Rep) %>% summarise(PercentDisrupt = mean(DistToClosestOncogene == -1, na.rm=T)) %>% .$PercentDisrupt %>% as.numeric()
    random_positions_median_distance[[gci]][[svi]] = this_random_positions %>% filter(!isDisruption) %>% group_by(Rep) %>% summarise(MedianDistance = median(DistToClosestOncogene)) %>% .$MedianDistance %>% as.numeric()
    
  }
  
}
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/PalmTreeTargetsCloseToOncogene.Rdata")



# ------------------------------------------------------------------
# -------------- ANALYSIS, STATISTICS, PLOTS -----------------------
# ------------------------------------------------------------------


rm(list=ls())
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/PalmTreeTargetsCloseToOncogene.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
svcallers = "Union"


# Annotate Tx with CircleGenome Classification (circle-circle, circle-genome, genome-genome)
txptinfo = txptinfo %>% 
  left_join(union_wgscircles %>% dplyr::select(BPID, CircleGenomeClass)) %>% 
  dplyr::rename(WGSCircles_CircleGenomeClass = CircleGenomeClass) %>% 
  left_join(union_circleseq %>% dplyr::select(BPID, CircleGenomeClass)) %>%
  dplyr::rename(CircleSeqCircles_CircleGenomeClass = CircleGenomeClass)


# Calculate statistics for All COSMIC genes vs. only Oncogenes vs. only TSGs
p_mediandist = list()
median_distance_circleseq_circlecircle = list()
for (gci in 1:length(geneclasses)){ 
  
  
  # For actual palm tree associated brekapoints, calculate distance to closest 
  # COSMIC gene and whether brekapoints disrupts a COSMIC gene
  txptinfo$DistToClosestOncogene = mapply(function (c, p) dist_to_closest_gene(c,p,geneclasses[[gci]]), txptinfo$TargetChrom, txptinfo$TargetPos)
  txptinfo$DisruptsCosmicOncogene = (txptinfo$DistToClosestOncogene == -1)
  
  
  # Go through the different Tx Callset, in this case only the merged 
  # Dataset (svcallers = "Union")
  p_mediandist[[gci]] = list()
  median_distance_circleseq_circlecircle[[gci]] = list()
  for (svi in 1:length(svcallers)){
    
    # Only one callset at a time
    this_txpt = txptinfo %>% filter(SVCaller == svcallers[svi])
    
    
    # Summary statistics for real breakpoints
    mean_distance[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>% distinct() %>% .$DistToClosestOncogene %>% mean()
    sd_distance[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>%  distinct() %>%.$DistToClosestOncogene %>% sd()
    median_distance[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>%  distinct() %>%.$DistToClosestOncogene %>% median()
    median_distance_circleseq_circlecircle[[gci]][[svi]] = 
      this_txpt %>% 
      filter(CircleSeqCircles_CircleGenomeClass == "circle-circle", 
             DistToClosestOncogene != -1, DistToClosestOncogene != Inf) %>%  distinct() %>%.$DistToClosestOncogene %>% median()
    
    
    # Summary statistics for randomly sampled breakpoints
    random_median_distances_ecdf = ecdf(random_positions_median_distance[[gci]][[svi]])
    p_mediandist[[gci]][[svi]] = random_median_distances_ecdf(median_distance[[gci]][[svi]])
    percent_disrupt[[gci]][[svi]] = mean(this_txpt$DistToClosestOncogene == -1)
    distance_distribution[[gci]][[svi]] = this_txpt %>% filter(DistToClosestOncogene != -1, DistToClosestOncogene != Inf)
    
    
    # Combine Real Data and Simulated Data into one Data Frame
    my_this_txpt = this_txpt %>% filter(DistToClosestOncogene != -1) %>%
      mutate(Chr = TargetChrom, Pos = TargetPos, Rep = 0, Draw = "Observed") %>%
      select(Chr, Pos, DistToClosestOncogene, Rep, Draw)
    my_random_positions = random_positions[[gci]][[svi]] %>%
      mutate(Pos = Start, Draw = "Random") %>%
      select(Chr, Pos, DistToClosestOncogene, Rep, Draw)
    my = rbind(my_this_txpt, my_random_positions)
    my$Draw = factor(my$Draw, levels=c("Random", "Observed"))
    
    
    # Beautify
    cols <- c("#000000", "#FF0000")
    cols.alpha <- c(grDevices::adjustcolor(cols[1], alpha.f = 0.01), cols[2])
    
    
    # Plot Distribution of Distances to closest COSMIC gene for real dataset vs. 
    # all randomized datasets
    my %>%
      filter(my$DistToClosestOncogene != -1) %>%
      arrange(Draw) %>%
      ggplot(aes(x=DistToClosestOncogene/1000, color=Draw, group=Rep)) +
      geom_density(size=0.5) + 
      scale_x_log10(labels = scales::comma) + 
      scale_color_manual(values = cols.alpha) + 
      theme_kons1() + 
      ggtitle(paste0(svcallers[svi],  "_", names(geneclasses)[gci])) +
      xlab("Distance to closest \n gene [kb]") + 
      ylab("Density") + 
      ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/DistToCosmic/", svcallers[svi], "_", names(geneclasses)[gci], "_DistanceDistribution.pdf"), height=3, width=3.5)
    
    
    # Plot Proportion of COSMIC gene-disrupting breakpoints for real vs. 
    # simulated data
    data.frame(x=random_positions_percent_disrupt[[gci]][[svi]]) %>%
      ggplot(aes(x=x)) + 
      geom_density() + 
      geom_vline(xintercept=percent_disrupt[[gci]][[svi]], color="red") +
      theme_kons1() + 
      ggtitle(paste0(svcallers[svi],  "_", names(geneclasses)[gci])) +
      xlab("Relative Proportion of \n Translocations that disrupt \n a gene \n Red = Observed Proportion") +
      ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/DistToCosmic/", svcallers[svi],  "_", names(geneclasses)[gci], "_Disruption.pdf"), height=3, width=3)
    
    
    # Same analysis but split by Circle-seq circle-genome classification of rearrangements
    my_this_txpt = this_txpt %>% filter(DistToClosestOncogene != -1) %>%
      mutate(Chr = TargetChrom, Pos = TargetPos, Rep = 0, Draw = CircleSeqCircles_CircleGenomeClass) %>%
      select(Chr, Pos, DistToClosestOncogene, Rep, Draw)
    my_this_txpt = my_this_txpt %>% filter(!is.na(Draw)) %>% mutate(Draw = droplevels(Draw))
    my_random_positions = random_positions[[gci]][[svi]] %>%
      mutate(Pos = Start, Draw = "Random") %>%
      select(Chr, Pos, DistToClosestOncogene, Rep, Draw)
    my = rbind(my_this_txpt, my_random_positions)
    my$Draw = factor(as.character(my$Draw), levels=c("Random", "circle-circle", "circle-genome", "genome-genome"))
    my %>%
      filter(DistToClosestOncogene != -1) %>%
      arrange(Draw) %>%
      #filter(Draw != "Random") %>% 
      ggplot(aes(x=DistToClosestOncogene/1000, color=Draw, group=interaction(Draw, Rep))) +
      geom_density(size=0.5) + 
      scale_x_log10(labels = scales::comma) + 
      scale_color_manual(values = c("Random" = grDevices::adjustcolor("black", alpha.f = 0.05),  "circle-circle"="firebrick2", "circle-genome"="forestgreen", "genome-genome"="steelblue2")) + 
      theme_kons1() + 
      ggtitle(paste0(svcallers[svi],  "_", names(geneclasses)[gci])) +
      xlab("Distance to closest \n gene [kb]") + 
      ylab("Density") + 
      ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/DistToCosmic/", svcallers[svi], "_", names(geneclasses)[gci], "_DistanceDistribution_CircleSeqCircleGenome.pdf"), height=3, width=3.5)
    circle_genome_classes = c("circle-circle", "circle-genome", "genome-genome")
    for (thisclass in circle_genome_classes){
      n_thisclass = this_txpt %>% filter(CircleSeqCircles_CircleGenomeClass == thisclass) %>% nrow()
      real_median_distance_thisclass = this_txpt %>% filter(WGSCircles_CircleGenomeClass == thisclass) %>% summarise(medianDist = median(DistToClosestOncogene)) %>% .$medianDist
      random_positions[[gci]][[svi]]$Index = rep(seq(1,max(as.numeric(random_positions[[gci]][[svi]]$Rep))), each=nrow(random_positions[[gci]][[svi]])/max(as.numeric(random_positions[[gci]][[svi]]$Rep)))
      random_median_distances_circlegenomeclass = random_positions[[gci]][[svi]] %>% filter(Index <= n_thisclass) %>% group_by(Rep) %>% summarise(medianDist = median(DistToClosestOncogene)) %>% .$medianDist
      random_median_distances_circlegenomeclass_ecdf = ecdf(random_median_distances_circlegenomeclass)
      p_mediandist_circlegenomeclass = random_median_distances_circlegenomeclass_ecdf(real_median_distance_thisclass)
      print(paste0(svcallers[svi], " ", names(geneclasses)[gci], " Circle-seq Circles ", thisclass, sprintf(" p=%.5f", p_mediandist_circlegenomeclass)))
    }
    
    
    # Same analysis but split by WGS circle-genome classification type of rearrangements
    my_this_txpt = this_txpt %>% filter(DistToClosestOncogene != -1) %>%
      mutate(Chr = TargetChrom, Pos = TargetPos, Rep = 0, Draw = WGSCircles_CircleGenomeClass) %>%
      select(Chr, Pos, DistToClosestOncogene, Rep, Draw)
    my_this_txpt = my_this_txpt %>% filter(!is.na(Draw)) %>% mutate(Draw = droplevels(Draw))
    my_random_positions = random_positions[[gci]][[svi]] %>%
      mutate(Pos = Start, Draw = "Random") %>%
      select(Chr, Pos, DistToClosestOncogene, Rep, Draw)
    my = rbind(my_this_txpt, my_random_positions)
    my$Draw = factor(as.character(my$Draw), levels=c("Random", "circle-circle", "circle-genome", "genome-genome"))
    my %>%
      filter(DistToClosestOncogene != -1) %>%
      ggplot(aes(x=DistToClosestOncogene/1000, color=Draw, group=interaction(Draw, Rep))) +
      geom_density(size=0.5, bw=0.5) + 
      scale_x_log10(labels = scales::comma) + 
      scale_color_manual(values = c("Random" = grDevices::adjustcolor("black", alpha.f = 0.05),  "circle-circle"="firebrick2", "circle-genome"="forestgreen", "genome-genome"="steelblue2")) + 
      theme_kons1() + 
      ggtitle(paste0(svcallers[svi],  "_", names(geneclasses)[gci])) +
      xlab("Distance to closest \n gene [kb]") + 
      ylab("Density") + 
      ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/DistToCosmic/", svcallers[svi], "_", names(geneclasses)[gci], "_DistanceDistribution_WGSCircleGenome.pdf"), height=3, width=3.5)
    circle_genome_classes = c("circle-circle", "circle-genome", "genome-genome")
    for (thisclass in circle_genome_classes){
      n_thisclass = this_txpt %>% filter(WGSCircles_CircleGenomeClass == thisclass) %>% nrow()
      real_median_distance_thisclass = this_txpt %>% filter(WGSCircles_CircleGenomeClass == thisclass) %>% summarise(medianDist = median(DistToClosestOncogene)) %>% .$medianDist
      random_positions[[gci]][[svi]]$Index = rep(seq(1,max(as.numeric(random_positions[[gci]][[svi]]$Rep))), each=nrow(random_positions[[gci]][[svi]])/max(as.numeric(random_positions[[gci]][[svi]]$Rep)))
      random_median_distances_circlegenomeclass = random_positions[[gci]][[svi]] %>% filter(Index <= n_thisclass) %>% group_by(Rep) %>% summarise(medianDist = median(DistToClosestOncogene)) %>% .$medianDist
      random_median_distances_circlegenomeclass_ecdf = ecdf(random_median_distances_circlegenomeclass)
      p_mediandist_circlegenomeclass = random_median_distances_circlegenomeclass_ecdf(real_median_distance_thisclass)
      print(paste0(svcallers[svi], " ", names(geneclasses)[gci], " WGS Circles ", thisclass, sprintf(" p=%.5f", p_mediandist_circlegenomeclass)))
    }
  }
}

