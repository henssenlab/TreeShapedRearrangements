rm(list=ls())

library(GenomicRanges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(biomaRt)

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/RNASeqData.Rdata")

txptinfo = txptinfo %>% 
  filter(PalmTreeChrom != TargetChrom) %>%
  filter(PalmTreeID != "NB2001_Novobreak_chr1:30634094-30685230") %>%
  as.data.frame()

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes <- c("hgnc_symbol", "start_position", "end_position", "gene_biotype")
filters <- c("chromosome_name","start","end")

txptinfo$Cohort = ifelse(grepl("NB2", txptinfo$Sample), "Berlin", "Peifer")

for (i in 1:nrow(txptinfo)){
  
  # Load Patient Metadata
  which_sample = txptinfo[i, "Sample"]
  print(which_sample)
  which_chr = txptinfo[i, "TargetChrom"] %>% gsub("chr", "", .)
  which_breakpoint = txptinfo[i, "TargetPos"]
  interval = 4000000 # size of interval of interest around breakpoint
  
  # use biomart to load all genes that are located in interval around breakpoint
  interval_start = which_breakpoint - floor(interval/2)
  interval_end = which_breakpoint + ceiling(interval/2)
  
  # only protein coding genes within that interval
  genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
    filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
    dplyr::select(hgnc_symbol, start_position, end_position)
  
  colnames(genes_within_interval) = c("Gene", "Start", "End")
  
  if (nrow(genes_within_interval) == 0) next
  
  means = inner_join(genes_within_interval, rna_tpm_means %>% filter(Cohort == txptinfo[i, "Cohort"]), by="Gene")
  tpmofinterest = inner_join(genes_within_interval, rna_tpm %>% filter(Cohort == txptinfo[i, "Cohort"], Sample == which_sample), by="Gene")
  
  if (nrow(tpmofinterest) == 0) next
  data = inner_join(means, tpmofinterest %>% spread(Sample, TPM))
  
  eval(parse(text=paste0("data$FoldChange = (data$", which_sample, "+ 0.00000001) / (data$MeanTPM+0.00000001)")))  
  data$LogFoldChange = log2(data$FoldChange)
  data$ModifiedZScore = ifelse((data$MadTPM == 0), 0*data[,which_sample], (data[,which_sample] - data$MedianTPM) / (data$MadTPM))
  data$Middle = (data$Start + data$End)/2
  
  # g=ggplot(data=data, aes(x=Middle, y=LogFoldChange, color=LogFoldChange)) + 
  #   geom_vline(xintercept=which_breakpoint, color="black", size=0.25) + 
  #   #geom_point(aes(x=Middle), size=1) + 
  #   geom_errorbarh(data=data, aes(xmin=Start, xmax=End), height=0.1) + 
  #   theme_minimal() + 
  #   geom_text_repel(aes(label=Gene), size=2) + 
  #   ggtitle(paste0("Palm Tree: ", gsub(":", "_", as.character(txptinfo[i, "PalmTreeID"])), "\nTarget: chr", as.character(which_chr), ":", as.character(which_breakpoint), "\nInterval around breakpoint: ", as.character(format(interval, scientific=FALSE)))) +
  #   xlab(paste0("Chr ", eval(which_chr))) + 
  #   ylab("TPM logFC") + 
  #   scale_y_continuous(limits=c(-10,5)) +
  #   scale_color_gradient2(midpoint=0, low="blue", mid="gray", high="red", na.value="red", limits=c(-10,5)) + 
  #   guides(color=FALSE)
  # gname = paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionPlots/EachTargetOnePlot/logFC/", gsub(":", "_", as.character(txptinfo[i, "PalmTreeID"])), "_", as.character(which_chr), "_", as.character(which_breakpoint), ".pdf")
  # ggsave(filename=gname, plot=g, width=12)
  
  g=ggplot(data=data, aes(x=(Middle-which_breakpoint)/1000, y=ModifiedZScore, color=ModifiedZScore)) + 
    geom_vline(xintercept=0, color="black", size=0.25, linetype="dashed") + 
    #geom_point(aes(x=Middle), size=1) + 
    geom_errorbarh(data=data, aes(xmin=(Start-which_breakpoint)/1000, xmax=(End-which_breakpoint)/1000), height=0.1) + 
    theme_classic() + 
    geom_text_repel(aes(label=Gene), size=2) + 
    ggtitle(paste0("Palm Tree: ", gsub(":", "_", as.character(txptinfo[i, "PalmTreeID"])), "\nTarget: chr", as.character(which_chr), ":", as.character(which_breakpoint))) +
    xlab("Distance from Breakpoint [kb]") + 
    ylab("Modified Z Score") + 
    geom_hline(yintercept = 0, linetype="dashed", size=0.25) + 
    #scale_y_continuous(limits=c(-10,5)) +
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


completedata = data.frame()
nodataadded = TRUE
for (i in 1:nrow(txptinfo)){
  
  # Load Patient Metadata
  which_sample = txptinfo[i, "Sample"]
  which_chr = txptinfo[i, "TargetChrom"] %>% gsub("chr", "", .)
  which_breakpoint = txptinfo[i, "TargetPos"]
  interval = 4000000 # size of interval of interest around breakpoint
  
  # use biomart to load all genes that are located in interval around breakpoint
  interval_start = which_breakpoint - floor(interval/2)
  interval_end = which_breakpoint + ceiling(interval/2)
  # only protein coding genes within that interval
  genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
    filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
    dplyr::select(hgnc_symbol, start_position, end_position)
  colnames(genes_within_interval) = c("Gene", "Start", "End")
  if (nrow(genes_within_interval) == 0) next
  
  means = inner_join(genes_within_interval, rna_tpm_means %>% filter(Cohort == txptinfo[i, "Cohort"]), by="Gene")
  tpmofinterest = inner_join(genes_within_interval, rna_tpm %>% filter(Cohort == txptinfo[i, "Cohort"], Sample == which_sample), by="Gene")
  
  if (nrow(tpmofinterest) == 0) next
  
  data = inner_join(means, tpmofinterest %>% spread(Sample, TPM))
  eval(parse(text=paste0("data$FoldChange = (data$", which_sample, "+ 0.000001) / (data$MeanTPM+0.000001)")))  
  eval(parse(text=paste0("data$ZScore = (data$", which_sample, "- data$MeanTPM) / (data$SdTPM+0.01)")))  
  data$Middle = (data$Start + data$End)/2
  data$ModifiedZScore = (data[,which_sample] - data$MedianTPM) / (data$MadTPM + 0.01)
  
  ngenesperside = 10
  leftdata = data %>% filter(Middle < which_breakpoint) %>% mutate(LeftRank = row_number(abs(Middle-which_breakpoint)))
  rightdata = data %>% filter(Middle >= which_breakpoint) %>% mutate(RightRank = row_number(abs(Middle-which_breakpoint)))
  data = full_join(leftdata, rightdata) %>% filter(LeftRank <= ngenesperside | RightRank <= ngenesperside)
  colnames(data)[colnames(data)==as.character(which_sample)] = "TPM"
  
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
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/PalmTreeXExpressionData.Rdata")
#load("~/Desktop/PalmTrees/Analysis/WorkspaceData/PalmTreeXExpressionData.Rdata")

# Export completedata to table
#completedata %>% arrange(desc(ModifiedZScore)) %>% View()
write.table(completedata %>% arrange(desc(ModifiedZScore)), "~/Desktop/PalmTrees/Results/Tables/CompleteExpressionAtInsertionSites.tsv", quote=F, col.names=T, row.names=F, sep="\t")

svcallers = "Union" #unique(txptinfo$SVCaller)
for (svi in 1:length(svcallers)){

  write.table(completedata %>% filter(SVCaller==svcallers[svi]) %>% arrange(desc(ModifiedZScore)), paste0("~/Desktop/PalmTrees/Results/Tables/", svcallers[svi], "_CompleteExpressionAtInsertionSites.tsv"), quote=F, col.names=T, row.names=F, sep="\t")

  
# Plot distributions of logFC and modified Z Score
completedata %>% 
  filter(SVCaller == svcallers[svi]) %>%
  filter(ModifiedZScore < 20, GeneDistFromBreakpoint<=10) %>%
  ggplot(aes(x=ModifiedZScore)) + geom_histogram(fill="steelblue", bins=50) + theme_classic() + geom_vline(xintercept=0, linetype="dashed") + xlab("Modified Z score") + ggtitle(svcallers[svi])
ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_DistributionOfZScores.pdf"), height=4, width=4)
# IMPORTANT: DISCUSS 0 VALUES (try bins=200 to see)
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

# Order Events (= Candidate Insertion Site = all Genes around one Target Breakpoint)

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


# FILTER OUT BREAKPOINTS THAT ARE WITHIN X OF ANOTHER BREAKPOINT TO AVOID DUPLICATE ROWS IN THE PLOT
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
svcallers = "Union" #unique(filtered_completedata$SVCaller)
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


svcallers = "Union"#unique(txptinfo$SVCaller)
for (svi in 1:length(svcallers)){
  txptinfo_this_svcaller = txptinfo %>% filter(SVCaller == svcallers[svi])
  permutated_logFC = list()
  actual_logFC = list()
  permutated_ModifiedZScore = list()
  actual_ModifiedZScore = list()
  
  for (i in 1:nrow(txptinfo_this_svcaller)){
    
    # Load Patient Metadata
    which_sample = txptinfo_this_svcaller[i, "Sample"]
    which_chr = txptinfo_this_svcaller[i, "TargetChrom"] %>% gsub("chr", "", .)
    which_breakpoint = txptinfo_this_svcaller[i, "TargetPos"]
    interval = 2000000 # size of interval of interest around breakpoint
    
    # use biomart to load all genes that are located in interval around breakpoint
    interval_start = which_breakpoint - floor(interval/2)
    interval_end = which_breakpoint + ceiling(interval/2)
    
    # only protein coding genes within that interval
    genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=as.character(which_chr),start=as.character(interval_start),end=as.character(interval_end)), mart=mart) %>% 
      filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
      dplyr::select(hgnc_symbol, start_position, end_position)
    colnames(genes_within_interval) = c("Gene", "Start", "End")
    
    if (nrow(genes_within_interval) == 0) next
    
    genes_within_interval$DistanceFromBreakpoint = abs(genes_within_interval$Start - which_breakpoint)
    genes_within_interval$DistanceFromBreakpointRank = rank(genes_within_interval$DistanceFromBreakpoint)
    genes_within_interval = genes_within_interval %>% filter(DistanceFromBreakpointRank<=5)
    
    means = inner_join(genes_within_interval, rna_tpm_means, by="Gene")
    tpmofinterest = inner_join(genes_within_interval, rna_tpm, by="Gene")
    
    data = inner_join(means %>% dplyr::select(Gene, Cohort, MeanTPM, SdTPM, MedianTPM, MadTPM), tpmofinterest, by="Gene")
    
    data$logFC = log2( (data$TPM + 0.000001) / (data$MeanTPM + 0.000001) )
    data$ModifiedZScore = (data$TPM - data$MedianTPM) / (data$MadTPM + 0.01)
    
    permutated_logFC[[i]] = data %>% filter(Sample != which_sample) %>% .$logFC
    actual_logFC[[i]] = data %>% filter(Sample == which_sample) %>% .$logFC
    
    permutated_ModifiedZScore[[i]] = data %>% filter(Sample != which_sample) %>% .$ModifiedZScore
    actual_ModifiedZScore[[i]] = data %>% filter(Sample == which_sample) %>% .$ModifiedZScore
    
    print(paste0(as.character(i), " / ", as.character(nrow(txptinfo_this_svcaller)), " of target sites are done."))
  }
  
  permutated_logFC = unlist(permutated_logFC)
  actual_logFC = unlist(actual_logFC)
  permutated_ModifiedZScore = unlist(permutated_ModifiedZScore)
  actual_ModifiedZScore = unlist(actual_ModifiedZScore)
  
  permoract =
    rbind(
      data.frame(isPerm = "AllOtherSamples", logFC = permutated_logFC, ModifiedZScore = permutated_ModifiedZScore),
      data.frame(isPerm = "RealSample", logFC = actual_logFC, ModifiedZScore = actual_ModifiedZScore)
    )
  
  ggplot(data=permoract, aes(x=logFC, color=isPerm, fill=isPerm)) + 
    geom_histogram(color=NA, bins=100) + 
    geom_vline(xintercept = 0, linetype="dashed") +
    facet_grid(isPerm~., scales = "free_y") +
    theme_classic() + 
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") 
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_logFC_PermVsActual.pdf"))
  
  ggplot(data=permoract, aes(x=logFC, color=isPerm)) + 
    geom_density() + 
    geom_vline(xintercept = 0, linetype="dashed") +
    theme_classic() + 
    scale_color_brewer(palette="Set1") 
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_logFC_PermVsActual_Density.pdf"))
  
  ggplot(data=permoract %>% filter(abs(ModifiedZScore) < 100), aes(x=ModifiedZScore, color=isPerm, fill=isPerm)) + 
    geom_histogram(bins=100, color=NA) + 
    geom_vline(xintercept = 0, linetype="dashed") +
    facet_grid(isPerm~., scales="free_y") +
    theme_classic() + 
    scale_color_brewer(palette="Set1") +
    scale_fill_brewer(palette="Set1") 
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_ModifiedZScore_PermVsActual.pdf"))
  
  ggplot(data=permoract %>% filter(abs(ModifiedZScore) < 100), aes(x=ModifiedZScore, color=isPerm)) + 
    geom_density() +
    geom_vline(xintercept = 0, linetype="dashed") +
    theme_classic() + 
    scale_color_brewer(palette="Set1")
  ggsave(paste0("~/Desktop/PalmTrees/Results/Figures/ExpressionTests/", svcallers[svi], "_ModifiedZScore_PermVsActual_Density.pdf"))
  
  save.image(paste0("~/Desktop/PalmTrees/Analysis/WorkspaceData/", svcallers[svi], "_ExpressionSampleRandomization.Rdata"))
  
}

# Make nicer plot

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



# ####################################################
# 
# ggplot(data=completedata%>%filter(abs(CoordDistFromBreakpoint)<5000000), aes(x=CoordDistFromBreakpoint, y=FoldChange)) + 
#   geom_point(alpha=0.1) + 
#   geom_vline(xintercept=0)+
#   geom_hline(yintercept=1)+
#   geom_density2d(color="steelblue") +
#   theme_minimal() + 
#   scale_y_log10() + 
#   theme(panel.grid = element_blank())
# 
# fc20data = completedata %>% filter(LeftRank > 25 | RightRank > 25) %>% group_by(EventXSide) %>% summarise(FoldChange20=mean(FoldChange))
# completedata_fc20sorted = full_join(completedata, fc20data)
# completedata_fc20sorted$EventXSide = reorder(completedata_fc20sorted$EventXSide, completedata_fc20sorted$FoldChange20, mean)
# ggplot(data=completedata_fc20sorted , aes(x=abs(TotalRank2), y=EventXSide)) + 
#   geom_tile(aes(fill=FoldChange)) + 
#   scale_fill_gradient2(midpoint=1, low="blue", mid="gray90", high="red", limits=c(0,20), na.value="red") +
#   theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
# 
# fc20data = completedata %>% filter(abs(TotalRank2) <= 25) %>% group_by(Event) %>% summarise(FoldChange20=mean(FoldChange))
# completedata_fc20sorted = full_join(completedata, fc20data)
# completedata_fc20sorted$Event = reorder(completedata_fc20sorted$Event, completedata_fc20sorted$FoldChange20, mean)
# ggplot(data=completedata_fc20sorted , aes(x=TotalRank2, y=Event)) + 
#   geom_tile(aes(fill=FoldChange)) + 
#   scale_fill_gradient2(midpoint=1, low="blue", mid="gray90", high="red", limits=c(0,20), na.value="red") +
#   theme(axis.text.y = element_blank(), axis.ticks.y=element_blank())
# 
# # sort by mean Fold change of the 50 closest genes to the breakpoint
# # then plot Events as rows with Gene Distance from Breakpoint as x an Fold Change as color
# fc20data = completedata %>% filter(abs(TotalRank2) <= 25) %>% group_by(Event) %>% summarise(FoldChange20=mean(FoldChange))
# completedata_fc20sorted = full_join(completedata, fc20data)
# completedata_fc20sorted$Event = reorder(completedata_fc20sorted$Event, completedata_fc20sorted$FoldChange20, mean)
# ggplot(data=completedata_fc20sorted , aes(x=CoordDistFromBreakpoint, y=Event)) + 
#   geom_point(aes(color=FoldChange), size=0.1) + 
#   scale_color_gradient2(midpoint=1, low="blue", mid="gray90", high="red", limits=c(0,10), na.value="red") +
#   theme(axis.text.y = element_blank(), axis.ticks.y=element_blank()) + 
#   geom_vline(xintercept = 0, size=0.2, linetype="dashed")
# 
# 
# 
# #############
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM)) %>% 
#   ggplot(aes(x=MeanTPM+0.000001)) + 
#   geom_histogram(bins=100) +
#   scale_x_log10()
# 
# rna_tpm %>% 
#   ggplot(aes(x=TPM+0.000001)) + 
#   geom_histogram(bins=100) +
#   scale_x_log10()
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=median(TPM)) %>% 
#   full_join(rna_tpm, by="Gene") %>%
#   ggplot(aes(x=(TPM+0.000001)/(MeanTPM+0.000001))) + 
#   geom_histogram(bins=100) + 
#   scale_x_log10()
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM)) %>% 
#   full_join(rna_tpm, by="Gene") %>%
#   ggplot(aes(x=(TPM+1)/(MeanTPM+1))) + 
#   geom_histogram(bins=100) + 
#   scale_x_log10()
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM)) %>% 
#   full_join(rna_tpm, by="Gene") %>%
#   ggplot(aes(x=(TPM+0.001)/(MeanTPM+0.001))) + 
#   geom_histogram(bins=100) + 
#   scale_x_log10()
# 
# # If I get rid of all TPM=0, then everything looks nice
# rna_tpm %>% 
#   filter(TPM>0) %>%
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM)) %>% 
#   full_join(rna_tpm %>% filter(TPM>0), by="Gene") %>%
#   ggplot(aes(x=(TPM+0.000001)/(MeanTPM+0.000001))) + 
#   geom_histogram(bins=100) + 
#   scale_x_log10()
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM+0.000001)) %>% 
#   full_join(rna_tpm, by="Gene") %>%
#   ggplot(aes(x=(TPM+0.000001)/(MeanTPM))) + 
#   geom_histogram(bins=100) + 
#   scale_x_log10()
# 
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM), SdTPM=sd(TPM)) %>% 
#   full_join(rna_tpm, by="Gene") %>%
#   ggplot(aes(x=(TPM-MeanTPM)/SdTPM)) + 
#   geom_histogram(bins=100) 
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM), SdTPM=sd(TPM)) %>% 
#   ggplot(aes(x=MeanTPM)) + 
#   geom_histogram()
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM), SdTPM=sd(TPM)) %>% 
#   ggplot(aes(x=SdTPM)) + 
#   geom_histogram()
# 
# rna_tpm %>% 
#   group_by(Gene) %>% 
#   summarise(MeanTPM=mean(TPM), SdTPM=sd(TPM)) %>% 
#   ggplot(aes(x=MeanTPM, y=SdTPM)) + 
#   geom_point() + 
#   scale_x_log10() + 
#   scale_y_log10()
# 

