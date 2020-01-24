rm(list=ls())
source("~/Desktop/PalmTrees/Analysis/Code/parseAllSV.R")
# system("gzip -d ~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/*/*/*.vcf.gz")

datapath = "~/Desktop/PalmTrees/Data/EliasMail190125/Konstantin_results_Variant_Callers_24.1.19/"
samples = dir(datapath)
samples = samples[!grepl("excluded", samples)]

delly = list()
smufin = list()
svaba = list()
for (i in 1:length(samples)){
  delly[[i]] = parse_delly_allsv(paste0(datapath, samples[i], "/delly2/delly2-svs.pre.vcf"), samples[i])
  smufin[[i]] = parse_smufin_allsv(paste0(datapath, samples[i], "/smufin/somatic_small_SVs.sample_", samples[i]), paste0(datapath, samples[i],"/smufin/somatic_large_SVs.sample_", samples[i]), samples[i])
  svaba[[i]] = parse_svaba_allsv(paste0(datapath, samples[i], "/svaba/svaba-indels.PASS.vcf"), paste0(datapath, samples[i],"/svaba/svaba-sv-fromindels.PASS.vcf"), paste0(datapath, samples[i],"/svaba/svaba-sv.PASS.vcf"), samples[i])
}
delly = do.call(rbind, delly) %>% filter(Quality == "PASS", PairedEndMedianMAPQ >= 20)
smufin = do.call(rbind, smufin)  %>% filter((Smufin.BNDCoverage >= 5) | ((Smufin.TinTumor >= 10) & (Smufin.TinControl <= 2)))

svaba = do.call(rbind, svaba) 
svaba_indel = svaba %>% filter(Class == "indel", Filter=="PASS") %>% as_tibble()
svaba_svfromindel = svaba %>% filter(Class == "indelfromsv", Filter=="PASS") %>% as_tibble()
svaba = svaba %>% filter(Class == "sv", Filter=="PASS") %>% as_tibble()
sv = bind_rows(delly, smufin, svaba)

smufin$Sample = gsub("CB", "NB", smufin$Sample)
delly$Sample = gsub("CB", "NB", delly$Sample)
svaba$Sample = gsub("CB", "NB", svaba$Sample)

sv = sv %>% filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
delly = delly %>% filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))
smufin = smufin %>% filter(!(ChrA == "chr2" & PosA >= 33000000 & PosA <= 34000000), !(ChrB == "chr2" & PosB >= 33000000 & PosB <= 34000000))

callinginfo = read.csv("~/Desktop/PalmTrees/Data/EliasMail190125/NB_samples_RUN_summary.csv", header=T, sep=";")
callinginfo$Sample = gsub("CB", "NB", callinginfo$Sample)

save.image("~/Desktop/PalmTrees/Analysis/WorkspaceData/AllSV.Rdata")


