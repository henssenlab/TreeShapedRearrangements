rm(list=ls())
library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)

# Load Palm Tree Data (run CallPalmTrees.R before)
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
txptinfo = txptinfo %>% as.data.frame()
txptinfo$Cohort = ifelse(grepl("NB2", txptinfo$Sample), "Berlin", "Peifer")


# Load parsed RNA-seq data (run ParseRNAseq.R before)
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/RNASeqData.Rdata")


# Get Gene Coordinates
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes <- c("hgnc_symbol", "start_position", "end_position", "gene_biotype")
filters <- c("chromosome_name","start","end")



# ------------------------------------------------------------------
# ----------- PLOT GENE EXPRESSION AROUND BREAKPOINTS --------------
# ------------------------------------------------------------------


# Go through all Palm Tree-associated Breakpoints (rows in txtptinfo)
for (i in 1:nrow(txptinfo)){
  
  which_sample = txptinfo[i, "Sample"]
  print(which_sample)
  
  
  # For each palm tree target site (i.e. palm tree branch), create an 
  # interval (+- 2Mb) around breakpoint to look for genes 
  which_chr = txptinfo[i, "TargetChrom"] %>% gsub("chr", "", .)
  which_breakpoint = txptinfo[i, "TargetPos"]
  interval = 4000000 
  interval_start = which_breakpoint - floor(interval/2)
  interval_end = which_breakpoint + ceiling(interval/2)
  
  
  # Consider only protein-coding genes within
  genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
    filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
    dplyr::select(hgnc_symbol, start_position, end_position)
  colnames(genes_within_interval) = c("Gene", "Start", "End")
  if (nrow(genes_within_interval) == 0) next
  
  
  # Get gene expression for genes around breakpoint (only matching cohort Berlin vs. Peifer)
  means = inner_join(genes_within_interval, rna_tpm_means %>% filter(Cohort == txptinfo[i, "Cohort"]), by="Gene")
  tpmofinterest = inner_join(genes_within_interval, rna_tpm %>% filter(Cohort == txptinfo[i, "Cohort"], Sample == which_sample), by="Gene")
  if (nrow(tpmofinterest) == 0) next
  
  
  # Calculate log Fold Change and Modified Z Score for genes around brekapoint 
  # with respect to all other samples in the cohort
  data = inner_join(means, tpmofinterest %>% spread(Sample, TPM))
  eval(parse(text=paste0("data$FoldChange = (data$", which_sample, "+ 0.00000001) / (data$MeanTPM+0.00000001)")))  
  data$LogFoldChange = log2(data$FoldChange)
  data$ModifiedZScore = ifelse((data$MadTPM == 0), 0*data[,which_sample], (data[,which_sample] - data$MedianTPM) / (data$MadTPM))
  data$Middle = (data$Start + data$End)/2
  
  
  # Plot all genes with respective expression around breakpoint
  g=ggplot(data=data, aes(x=(Middle-which_breakpoint)/1000, y=ModifiedZScore, color=ModifiedZScore)) + 
    geom_vline(xintercept=0, color="black", size=0.25, linetype="dashed") + 
    geom_errorbarh(data=data, aes(xmin=(Start-which_breakpoint)/1000, xmax=(End-which_breakpoint)/1000), height=0.1) + 
    theme_classic() + 
    geom_text_repel(aes(label=Gene), size=2) + 
    ggtitle(paste0("Palm Tree: ", gsub(":", "_", as.character(txptinfo[i, "PalmTreeID"])), "\nTarget: chr", as.character(which_chr), ":", as.character(which_breakpoint))) +
    xlab("Distance from Breakpoint [kb]") + 
    ylab("Modified Z Score") + 
    geom_hline(yintercept = 0, linetype="dashed", size=0.25) + 
    scale_x_continuous(limits=c(-1200, 1200)) + 
    scale_color_gradient2(midpoint=0, low="red", mid="gray", high="red", na.value="red", limits=c(-5,5)) + 
    guides(color=FALSE) +
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5))
  gname = paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionPlots/EachTargetOnePlot/ModifiedZ/", txptinfo[i,"SVCaller"], "/", gsub(":", "_", as.character(txptinfo[i, "PalmTreeID"])), "_", as.character(which_chr), "_", as.character(which_breakpoint), ".pdf")
  ggsave(filename=gname, plot=g, width=12)
  
}



# ------------------------------------------------------------------
# ---------------------- EXPRESSION HEATMAP ------------------------
# ------------------------------------------------------------------


completedata = data.frame()
nodataadded = TRUE
for (i in 1:nrow(txptinfo)){
  
  which_sample = txptinfo[i, "Sample"]
  
  
  # For each palm tree target site (i.e. palm tree branch), create an 
  # interval (+- 2Mb) around breakpoint to look for genes 
  which_chr = txptinfo[i, "TargetChrom"] %>% gsub("chr", "", .)
  which_breakpoint = txptinfo[i, "TargetPos"]
  interval = 4000000 # size of interval of interest around breakpoint
  interval_start = which_breakpoint - floor(interval/2)
  interval_end = which_breakpoint + ceiling(interval/2)
  
  
  # Get protein-coding genes within interval
  genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
    filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
    dplyr::select(hgnc_symbol, start_position, end_position)
  colnames(genes_within_interval) = c("Gene", "Start", "End")
  if (nrow(genes_within_interval) == 0) next
  
  
  # Get gene expression for genes around breakpoint (only matching cohort Berlin vs. Peifer)
  means = inner_join(genes_within_interval, rna_tpm_means %>% filter(Cohort == txptinfo[i, "Cohort"]), by="Gene")
  tpmofinterest = inner_join(genes_within_interval, rna_tpm %>% filter(Cohort == txptinfo[i, "Cohort"], Sample == which_sample), by="Gene")
  if (nrow(tpmofinterest) == 0) next
  
  
  # Calculate log Fold Change and Modified Z Score for genes around brekapoint 
  # with respect to all other samples in the cohort
  data = inner_join(means, tpmofinterest %>% spread(Sample, TPM))
  eval(parse(text=paste0("data$FoldChange = (data$", which_sample, "+ 0.000001) / (data$MeanTPM+0.000001)")))  
  eval(parse(text=paste0("data$ZScore = (data$", which_sample, "- data$MeanTPM) / (data$SdTPM+0.01)")))  
  data$Middle = (data$Start + data$End)/2
  data$ModifiedZScore = (data[,which_sample] - data$MedianTPM) / (data$MadTPM + 0.01)
  
  
  # Consider only the 10 closest genes to breakpoint
  ngenesperside = 10
  leftdata = data %>% filter(Middle < which_breakpoint) %>% mutate(LeftRank = row_number(abs(Middle-which_breakpoint)))
  rightdata = data %>% filter(Middle >= which_breakpoint) %>% mutate(RightRank = row_number(abs(Middle-which_breakpoint)))
  data = full_join(leftdata, rightdata) %>% filter(LeftRank <= ngenesperside | RightRank <= ngenesperside)
  colnames(data)[colnames(data)==as.character(which_sample)] = "TPM"
  
  
  # Create several columns that will be used for plotting later
  data$TotalRank = mapply(function (lr, rr) ifelse(is.na(rr), -lr, rr), data$LeftRank, data$RightRank)
  data$GeneDistFromBreakpoint = mapply(function (lr, rr) ifelse(is.na(rr), lr, rr), data$LeftRank, data$RightRank)
  data$Sample = as.character(which_sample)
  data$Chromosome = which_chr
  data$Breakpoint = which_breakpoint
  data$CoordDistFromBreakpoint = data$Breakpoint-data$Middle
  data$Event = interaction(data$Sample, data$Chromosome, data$Breakpoint) %>% droplevels()
  data$Side = mapply(function (lr, rr) ifelse(is.na(rr), "LeftSide", "RightSide"), data$LeftRank, data$RightRank)
  data$EventXSide = interaction(data$Event, data$Side) %>% droplevels()
  data$SVCaller = as.character(txptinfo[i, "SVCaller"])
  data$PalmTreeID = as.character(txptinfo[i, "PalmTreeID"])
  
  
  if (nodataadded){
    completedata = data
    nodataadded = FALSE
  } else {
    completedata = rbind(completedata, data)
  }
  
}


# Save completedata
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/PalmTreeXExpressionData.Rdata")
write.table(completedata %>% arrange(desc(ModifiedZScore)), "~/Desktop/PalmTrees/Results/Tables/CompleteExpressionAtInsertionSites.tsv", quote=F, col.names=T, row.names=F, sep="\t")


# Loop through callsets separately (we only consider the merged callset)
svcallers = "Union"
for (svi in 1:length(svcallers)){
  
  
  # Write expression data for this call set
  write.table(completedata %>% filter(SVCaller==svcallers[svi]) %>% arrange(desc(ModifiedZScore)), paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_CompleteExpressionAtInsertionSites.tsv"), quote=F, col.names=T, row.names=F, sep="\t")

  
  # Plot distributions of Z score
  completedata %>% 
    filter(SVCaller == svcallers[svi]) %>%
    filter(ModifiedZScore < 20, GeneDistFromBreakpoint<=10) %>%
    ggplot(aes(x=ModifiedZScore)) + geom_histogram(fill="steelblue", bins=50) + theme_classic() + geom_vline(xintercept=0, linetype="dashed") + xlab("Modified Z score") + ggtitle(svcallers[svi])
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_DistributionOfZScores.pdf"), height=4, width=4)

  
  # Plot distributions of log fold change
  completedata %>% 
    filter(SVCaller ==   svcallers[svi]) %>%
    ggplot(aes(x=log2(FoldChange))) + geom_histogram(fill="steelblue", bins=50) + theme_classic() + geom_vline(xintercept=0, linetype="dashed") + xlab("Log2 Fold Change") + ggtitle(svcallers[svi])
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_DistributionLogFC.pdf"), height=4, width=4)
  
  
  # Recurrently affected genes around breakpoints
  completedata %>% 
    filter(SVCaller ==  svcallers[svi]) %>%
    dplyr::select(Gene, Sample, Chromosome, Start, End) %>% distinct() %>% 
    group_by(Gene, Chromosome, Start, End) %>% summarise(n=n()) %>% 
    arrange(desc(n)) %>% 
    write.table(paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_GenesWithin1MbOfInsertion_Recurrence.tsv"), quote=F, col.names=T, row.names=F, sep="\t")
  
  
  # Create heatmap for gene expression around all palm tree branch sites
  completedata %>% 
    filter(SVCaller ==  svcallers[svi]) %>%
    full_join(completedata %>% filter(SVCaller==svcallers[svi], GeneDistFromBreakpoint <= 10) %>% group_by(Event) %>% summarise(MeanModZ=mean(ModifiedZScore)), by="Event") %>%
    mutate(Event = factor(completedata %>% filter(SVCaller==svcallers[svi]) %>% .$Event, levels = completedata %>% filter(SVCaller==svcallers[svi]) %>% arrange(MeanModZ) %>% .$Event %>% unique())) %>%
    filter(GeneDistFromBreakpoint <= 10) %>%
    ggplot(aes(x=TotalRank, y=Event)) + 
    geom_vline(xintercept=0, size=0.25, linetype="dashed") + 
    geom_tile(aes(fill=ModifiedZScore)) + 
    scale_fill_gradient2(name="Modified Z Score", midpoint=0, low="blue", mid="gray", high="red", limits=c(-4,4), na.value="red") +
    theme(axis.text = element_blank(), axis.ticks=element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + 
    xlab("+/- 10 Genes around Breakpoint") +
    ylab("Candidate Integration Sites") +
    ggtitle(paste0("Expression Change around Candidate Reintegration Sites\n", svcallers[svi])) +
    theme(text = element_text(family = "Helvetica"),
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5),
          legend.position="bottom")
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_ExpressionHeatmap.pdf"), height=5, width=5)
}


# Filter out breakpoints that are very close to each other which 
# would produce a few duplicated lines in the heatmap
completedata$PalmTreeXTargetChr = paste0(completedata$PalmTreeID, ".Target.chr", completedata$Chromosome)
palmtreeXtargetchr = unique(completedata$PalmTreeXTargetChr)
filtered_completedata = list()
for (i in 1:length(palmtreeXtargetchr)){
  this = completedata %>% filter(PalmTreeXTargetChr == palmtreeXtargetchr[i]) 
  this_breakpoints = unique(this$Breakpoint)
  selected_breakpoints = this_breakpoints[!duplicated(ceiling(this_breakpoints / 4000000))]
  filtered_completedata[[i]] = this %>% filter(Breakpoint %in% selected_breakpoints)
}
filtered_completedata = do.call(rbind, filtered_completedata)

# Number of tumors and breakpoints in heatmap to put in figure legend
filtered_completedata %>% 
  filter(SVCaller == "Union") %>%
  full_join(filtered_completedata %>% filter(SVCaller=="Union", GeneDistFromBreakpoint <= 10) %>% group_by(Event) %>% summarise(MeanModZ=mean(ModifiedZScore)), by="Event") %>%
  mutate(Event = factor(filtered_completedata %>% filter(SVCaller=="Union") %>% .$Event, levels = filtered_completedata %>% filter(SVCaller=="Union") %>% arrange(MeanModZ) %>% .$Event %>% unique())) %>%
  filter(GeneDistFromBreakpoint <= 10) %>%
  summarise(nEvents = n_distinct(Event),
            nSamples = n_distinct(Sample))

# Redo the heatmap, same as above
svcallers = "Union"
for (svi in 1:length(svcallers)){
  filtered_completedata %>% 
    filter(SVCaller ==  svcallers[svi]) %>%
    full_join(filtered_completedata %>% filter(SVCaller==svcallers[svi], GeneDistFromBreakpoint <= 10) %>% group_by(Event) %>% summarise(MeanModZ=mean(ModifiedZScore)), by="Event") %>%
    mutate(Event = factor(filtered_completedata %>% filter(SVCaller==svcallers[svi]) %>% .$Event, levels = filtered_completedata %>% filter(SVCaller==svcallers[svi]) %>% arrange(MeanModZ) %>% .$Event %>% unique())) %>%
    filter(GeneDistFromBreakpoint <= 10) %>%
    ggplot(aes(x=TotalRank, y=Event)) + 
    geom_vline(xintercept=0, size=0.25, linetype="dashed") + 
    geom_tile(aes(fill=ModifiedZScore)) + 
    scale_fill_gradient2(name="Modified Z Score", midpoint=0, low="blue", mid="gray", high="red", limits=c(-4,4), na.value="red") +
    theme(axis.text = element_blank(), axis.ticks=element_blank(), panel.grid = element_blank(), panel.background = element_blank()) + 
    xlab("+/- 10 Genes around Breakpoint") +
    ylab("Candidate Integration Sites") +
    ggtitle(paste0("Expression Change around Candidate Reintegration Sites\n", svcallers[svi])) +
    theme(text = element_text(family = "Helvetica"),
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 12, hjust=0.5),
          legend.position="bottom")
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_ExpressionHeatmapFiltered.pdf"), height=5, width=5)
}



# ------------------------------------------------------------------
# -------- COMPARE ACTUAL EXPRESSION AROUND BREAKPOINTS VS. --------
# ------------ PERMUTING EXPRESSION VALUES OVER SAMPLES ------------
# ------------------------------------------------------------------

svcallers = "Union"
for (svi in 1:length(svcallers)){
  
  txptinfo_this_svcaller = txptinfo %>% filter(SVCaller == svcallers[svi])
  permutated_logFC = list()
  actual_logFC = list()
  permutated_ModifiedZScore = list()
  actual_ModifiedZScore = list()
  
  for (i in 1:nrow(txptinfo_this_svcaller)){
    
    which_sample = txptinfo_this_svcaller[i, "Sample"]
   
    
    # For each palm tree target site (i.e. palm tree branch), create an 
    # interval (+- 2Mb) around breakpoint to look for genes 
    which_chr = txptinfo_this_svcaller[i, "TargetChrom"] %>% gsub("chr", "", .)
    which_breakpoint = txptinfo_this_svcaller[i, "TargetPos"]
    interval = 2000000 # size of interval of interest around breakpoint
    interval_start = which_breakpoint - floor(interval/2)
    interval_end = which_breakpoint + ceiling(interval/2)
    
    
    # Get protein-coding genes within interval
    genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
      filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
      dplyr::select(hgnc_symbol, start_position, end_position)
    colnames(genes_within_interval) = c("Gene", "Start", "End")
    if (nrow(genes_within_interval) == 0) next
    
    
    # Compute distance from breakpoint for all genes, only consider 5 closest genes 
    genes_within_interval$DistanceFromBreakpoint = abs(genes_within_interval$Start - which_breakpoint)
    genes_within_interval$DistanceFromBreakpointRank = rank(genes_within_interval$DistanceFromBreakpoint)
    genes_within_interval = genes_within_interval %>% filter(DistanceFromBreakpointRank<=5)
    
    
    # Calculate log Fold Change and Modified Z Score for genes around brekapoint 
    means = inner_join(genes_within_interval, rna_tpm_means, by="Gene")
    tpmofinterest = inner_join(genes_within_interval, rna_tpm, by="Gene")
    data = inner_join(means %>% dplyr::select(Gene, Cohort, MeanTPM, SdTPM, MedianTPM, MadTPM), tpmofinterest, by="Gene")
    data$logFC = log2( (data$TPM + 0.000001) / (data$MeanTPM + 0.000001))
    data$ModifiedZScore = (data$TPM - data$MedianTPM) / (data$MadTPM + 0.01)
    
    
    # Save log FC and modified Z score of the genes of interest for all other samples.
    # This will allow to permutate the expression values to assess significance of 
    # expression change close to breakpoints
    permutated_logFC[[i]] = data %>% filter(Sample != which_sample) %>% .$logFC
    actual_logFC[[i]] = data %>% filter(Sample == which_sample) %>% .$logFC
    permutated_ModifiedZScore[[i]] = data %>% filter(Sample != which_sample) %>% .$ModifiedZScore
    actual_ModifiedZScore[[i]] = data %>% filter(Sample == which_sample) %>% .$ModifiedZScore
    
    print(paste0(as.character(i), " / ", as.character(nrow(txptinfo_this_svcaller)), " of target sites are done."))
  
  }
  
  
  # Combine permuted expression datasets and the actual dataset into one table
  # for convenient plotting
  permutated_logFC = unlist(permutated_logFC)
  actual_logFC = unlist(actual_logFC)
  permutated_ModifiedZScore = unlist(permutated_ModifiedZScore)
  actual_ModifiedZScore = unlist(actual_ModifiedZScore)
  permoract =
    rbind(
      data.frame(isPerm = "AllOtherSamples", logFC = permutated_logFC, ModifiedZScore = permutated_ModifiedZScore),
      data.frame(isPerm = "RealSample", logFC = actual_logFC, ModifiedZScore = actual_ModifiedZScore)
    )
  
  # Plot histogram of fold change of genes around a breakpoint for the actual
  # versus for the expression-permuted datasets
  ggplot(data=permoract, aes(x=logFC, color=isPerm, fill=isPerm)) + 
    geom_histogram(color=NA, bins=100) + 
    geom_vline(xintercept = 0, linetype="dashed") +
    facet_grid(isPerm~., scales = "free_y") +
    theme_classic() + 
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") 
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_logFC_PermVsActual.pdf"))
  
  
  # Plot density
  ggplot(data=permoract, aes(x=logFC, color=isPerm)) + 
    geom_density() + 
    geom_vline(xintercept = 0, linetype="dashed") +
    theme_classic() + 
    scale_color_brewer(palette="Set1") 
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_logFC_PermVsActual_Density.pdf"))
  
  
  # Plot histogram of modified Z scores of genes around a breakpoint for the actual
  # versus for the expression-permuted datasets
  ggplot(data=permoract %>% filter(abs(ModifiedZScore) < 100), aes(x=ModifiedZScore, color=isPerm, fill=isPerm)) + 
    geom_histogram(bins=100, color=NA) + 
    geom_vline(xintercept = 0, linetype="dashed") +
    facet_grid(isPerm~., scales="free_y") +
    theme_classic() + 
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") 
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_ModifiedZScore_PermVsActual.pdf"))
  
  
  # Plot density
  ggplot(data=permoract %>% filter(abs(ModifiedZScore) < 100), aes(x=ModifiedZScore, color=isPerm)) + 
    geom_density() +
    geom_vline(xintercept = 0, linetype="dashed") +
    theme_classic() + 
    scale_color_brewer(palette="Set1")
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_ModifiedZScore_PermVsActual_Density.pdf"))
  
  
  # Save permutation data
  save.image(paste0("~/Desktop/PalmTrees/Analysis/WorkspaceData/", svcallers[svi], "_ExpressionSampleRandomization.Rdata"))
}



# ------------------------------------------------------------------
# -------------- PRETTY PLOTTING FOR MANUSCRIPT --------------------
# ------------------------------------------------------------------

# Make small bar plots for the manuscript
completedata$EventXCaller = paste0(completedata$Event, ".", completedata$SVCaller)
events = unique(completedata$EventXCaller)
for (i in 1:length(events)){
  this_insertion = completedata %>% filter(EventXCaller == events[i]) %>% filter(abs(TotalRank)<=5)
  this_insertion %>%
    dplyr::select(TotalRank, ModifiedZScore, Gene) %>% 
    distinct() %>% 
    ggplot(aes(x=TotalRank, y=ModifiedZScore, fill=ModifiedZScore)) + 
    geom_col() + 
    geom_vline(xintercept = 0, linetype="dashed", size=0.25) + 
    theme_classic() +
    scale_x_continuous(breaks = this_insertion$TotalRank,
                       labels = this_insertion$Gene) + 
    theme(text = element_text(family = "Helvetica"),
          axis.text.y = element_text(color="black"),
          axis.text.x = element_text(face = "italic", color = "black", angle=45,vjust = 1, hjust = 1),
          axis.text = element_text(size = 10), 
          axis.title =  element_text(size = 10), 
          plot.title = element_text(face="plain", size = 6, hjust=0.5)) +
    xlab("") +
    ylab("Modified Z Score") +
    ggtitle(paste0("Palm Tree: ", gsub(":", "_", as.character(this_insertion[1, "PalmTreeID"])), "\nTarget: chr", as.character(this_insertion[1, "Chromosome"]), ":", as.character(this_insertion[1,"Breakpoint"]))) +
    guides(fill=F)+
    scale_fill_gradient2(name="Modified Z Score", midpoint=0, low="blue", mid="gray95", high="red", limits=c(min(c(-3,this_insertion$ModifiedZScore)),max(c(3,this_insertion$ModifiedZScore))), na.value="red") +
    ggsave(filename=paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionPlots/EachTargetOnePlot/ModifiedZColPlots/", this_insertion[1,"SVCaller"], "/", gsub(":", "_", as.character(this_insertion[1, "PalmTreeID"])), "_", as.character(this_insertion[1, "Chromosome"]), "_", as.character(this_insertion[1, "Breakpoint"]), ".pdf"), height=3, width=3)
}
