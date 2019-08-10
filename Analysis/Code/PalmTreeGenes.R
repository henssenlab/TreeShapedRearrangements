library(GenomicRanges)
library(IRanges)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyr)
library(biomaRt)


# Load Palm Tree Data
load("~/Desktop/PalmTrees/Analysis/WorkspaceData/saved_callPalmTree.Rdata")


# Load Gene Coordinates
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")
filters <- c("chromosome_name","start","end")


# Iterate over palm tree regions
for (i in 1:nrow(palmtrees)){
  
  
  # Get all genes within current palm tree region
  genes_within_interval = getBM(attributes=attributes, filters=filters, values=list(chromosome=palmtrees[[i,"Chr"]]%>% gsub("chr", "", .),start=as.character(palmtrees[[i,"FirstElement"]]),end=as.character(palmtrees[[i,"LastElement"]])), mart=mart) %>% filter(hgnc_symbol != "", gene_biotype == "protein_coding")
  colnames(genes_within_interval) = c("Gene", "Chr", "GeneStart", "GeneEnd", "GeneBiotype")
  if (nrow(genes_within_interval)==0) next
  genes_within_interval$Chr = paste0("chr", as.character(genes_within_interval$Chr))
  genes_within_interval$PalmTreeID = palmtrees[[i,"PalmTreeID"]]
  genes_within_interval$Sample = palmtrees[[i,"Sample"]]
  genes_within_interval$SVCaller = palmtrees[[i,"SVCaller"]]
  if (i==1){
    all_genes_within_palmtrees = genes_within_interval
  }else{
    all_genes_within_palmtrees = bind_rows(all_genes_within_palmtrees, genes_within_interval)
  }

}


# All Genes with Sample information
write.table(x=all_genes_within_palmtrees, file = "~/Desktop/PalmTrees/Results/Tables/GenesInPalmTrees/GenesInPalmTreesBySample.tsv", quote=F, sep = "\t", row.names=F)


# Compute which genes are recurrently located in palm tree regions  
all_genes_within_palmtrees %>%
  dplyr::select(Gene, Sample, SVCaller, Chr, GeneStart, GeneEnd) %>%
  distinct() %>% group_by(Gene, SVCaller) %>%
  summarise(NumberOfPalmTrees=n(), Chr=Chr[[1]], GeneStart=GeneStart[[1]], GeneEnd=GeneEnd[[1]]) %>% arrange(desc(NumberOfPalmTrees)) %>%
  write.table(x=., file = "~/Desktop/PalmTrees/Results/Tables/GenesInPalmTrees/GeneRecurrenceInPalmTree.tsv", quote=F, sep = "\t", row.names=F)


# Compute which genes are recurrently located in palm tree regions
# Analyis split up by Tx call set
svcallers = unique(palmtrees$SVCaller)
for (i in 1:length(svcallers)){
  all_genes_within_palmtrees %>%
    filter(SVCaller == svcallers[i]) %>%
    dplyr::select(Gene, Sample, SVCaller, Chr, GeneStart, GeneEnd) %>%
    distinct() %>% group_by(Gene, SVCaller) %>%
    summarise(NumberOfPalmTrees=n(), Chr=Chr[[1]], GeneStart=GeneStart[[1]], GeneEnd=GeneEnd[[1]]) %>% arrange(desc(NumberOfPalmTrees)) %>%
    write.table(x=., file = paste0("~/Desktop/PalmTrees/Results/Tables/GenesInPalmTrees/", svcallers[i], "_GeneRecurrenceInPalmTree.tsv"), quote=F, sep = "\t", row.names=F)
}


# Import list of COSMIC genes
cosmic = read.table("~/Desktop/PalmTrees/Data/COSMIC_Census_allTue Dec 11 07_52_32 2018.tsv", header=T, sep="\t")
all_genes_within_palmtrees %>%
  filter(Gene %in% cosmic$Gene.Symbol) %>%
  write.table(x=., file = "~/Desktop/PalmTrees/Results/Tables/GenesInPalmTrees/GenesInPalmTreesBySample_OnlyCosmic.tsv", quote=F, sep = "\t", row.names=F)


# Compute which COSMIC genes are recurrently located in palm tree regions  
all_genes_within_palmtrees %>%
  dplyr::select(Gene, Sample, SVCaller, Chr, GeneStart, GeneEnd) %>%
  filter(Gene %in% cosmic$Gene.Symbol) %>%
  distinct() %>% group_by(Gene, SVCaller) %>%
  summarise(NumberOfPalmTrees=n(), Chr=Chr[[1]], GeneStart=GeneStart[[1]], GeneEnd=GeneEnd[[1]]) %>% arrange(desc(NumberOfPalmTrees)) %>%
  write.table(x=., file = "~/Desktop/PalmTrees/Results/Tables/GenesInPalmTrees/GeneRecurrenceInPalmTree_OnlyCosmic.tsv", quote=F, sep = "\t", row.names=F)


# Compute which COSMIC genes are recurrently located in palm tree regions
# Analyis split up by Tx call set
svcallers = unique(palmtrees$SVCaller)
for (i in 1:length(svcallers)){
  all_genes_within_palmtrees %>%
    filter(SVCaller == svcallers[i]) %>%
    filter(Gene %in% cosmic$Gene.Symbol) %>%
    dplyr::select(Gene, Sample, SVCaller, Chr, GeneStart, GeneEnd) %>%
    distinct() %>% group_by(Gene, SVCaller) %>%
    summarise(NumberOfPalmTrees=n(), Chr=Chr[[1]], GeneStart=GeneStart[[1]], GeneEnd=GeneEnd[[1]]) %>% arrange(desc(NumberOfPalmTrees)) %>%
    write.table(x=., file = paste0("~/Desktop/PalmTrees/Results/Tables/GenesInPalmTrees/", svcallers[i], "_GeneRecurrenceInPalmTree_OnlyCosmic.tsv"), quote=F, sep = "\t", row.names=F)
}

