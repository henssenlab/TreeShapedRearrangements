rm(list=ls())
library(dplyr)
library(biomaRt)
library(parallel)
cores = detectCores()
source("~/Desktop/PalmTrees/Analysis/Code/hasOverlap.R")
source("~/Desktop/PalmTrees/Analysis/Code/ParseSVCallerData.R")

load("~/Desktop/PalmTrees/Analysis/WorkspaceData/ascat_cnv.Rdata")
rm(list=c("cnv_gr", "cnv"))
ascat_cnv$isLOH = (ascat_cnv$TumorMinorAlleleCopyNumber == 0) & (ascat_cnv$NormalMinorAlleleCopyNumber == 1)

mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
attributes <- c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "gene_biotype")

genes = 
  getBM(attributes=attributes, mart=mart) %>% 
  filter(hgnc_symbol != "", gene_biotype == "protein_coding") %>%
  mutate(Gene=hgnc_symbol, GeneChr = paste0("chr", chromosome_name), GeneStart = start_position, GeneEnd = end_position) %>%
  dplyr::select(Gene, GeneChr, GeneStart, GeneEnd) %>% 
  filter(isdefchrom2(GeneChr)) %>%
  as_tibble()

samples = unique(ascat_cnv$Sample)
gene_cn = mclapply(1:length(samples), function(i){
  this_genes = genes %>% mutate(Sample = samples[i], CopyNumber = 2, LOH=F)
  this_cnv = ascat_cnv %>% filter(Sample == samples[i])
  for (j in 1:nrow(this_cnv)){
    this_genes[hasOverlap_withChr(this_cnv[[j, "Chr"]], this_cnv[[j, "Start"]], this_cnv[[j, "End"]], this_genes$GeneChr, this_genes$GeneStart, this_genes$GeneEnd),
               "CopyNumber"] = this_cnv[[j, "TumorTotalCopyNumber"]]
    this_genes[hasOverlap_withChr(this_cnv[[j, "Chr"]], this_cnv[[j, "Start"]], this_cnv[[j, "End"]], this_genes$GeneChr, this_genes$GeneStart, this_genes$GeneEnd),
               "CopyNumberFC"] = this_cnv[[j, "TumorTotalCopyNumber"]] / this_cnv[[j, "NormalTotalCopyNumber"]]
    this_genes[hasOverlap_withChr(this_cnv[[j, "Chr"]], this_cnv[[j, "Start"]], this_cnv[[j, "End"]], this_genes$GeneChr, this_genes$GeneStart, this_genes$GeneEnd),
               "LOH"] = this_cnv[[j, "isLOH"]]
  }
  return(this_genes)
},
mc.cores=cores)
gene_cn = do.call(rbind, gene_cn)
rm(list=setdiff(ls(), c("ascat_cnv", "gene_cn")))
save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/GeneCopyNumber.Rdata")
