This repository contains code and data from Koche et al. [Extrachromosomal circular DNA drives oncogenic genome remodeling in neuroblastoma](https://doi.org/10.1038/s41588-019-0547-z) *Nat Genet* (2020).

### Data

**Circle Calls**
- `Data/CircleseqCircles.bed` contains Circle-seq circle calls for our dataset
- `Data/WGSinferredCircles_eccDNA.txt` contains WGS-inferred circle calls for our dataset (eccDNA only)
- `Data/WGSinferredCircles_ecDNA.txt` contains WGS-inferred circle calls for our dataset (ecDNA only)

**WGS Rearrangement Calls**
- `Data/MergedSV_StrictestFiltering_CircleSeqCircleAnnotation.txt` are the merged and filtered interchromosomal rearrangement calls (see Structural Variant Merging) which were overlapped with Circle-seq circle calls for classification *(available upon request)*
- `Data/MergedSV_StrictestFiltering_WGSCircleAnnotation.txt` are the merged and filtered interchromosomal rearrangement calls (see Structural Variant Merging) which were overlapped with WGS-inferred circle calls for classification *(available upon request)*
- `Data/PedPanCanSVs.csv` contains structural variants in the DKFZ Pediatric Pan-Cancer dataset published by Gröbner et al. 2018 (data downloaded from the corresponding [R2 platform](https://hgserver1.amc.nl/cgi-bin/r2/main.cgi?&dscope=DKFZ_PED&option=about_dscope))
- `Data/PedPanCanMeta.csv` contains metadata such as the tumor entity for the DKFZ Pediatric Pan-Cancer dataset

**Expression Data**
- `Data/berlin_cohort_rnaseq_fpkm.txt` is RNA-seq data for the Berlin neuroblastoma cohort *(available upon request)*
- `Data/peifer_54nb_fpkms.txt` is RNA-seq data published by Peifer et al.

**Clinical Data**
- `Data/ClinicalData.csv` specifies risk group and survival data

**Other**
- `Data/BroadIntervals/b37-[...]_MinusBlckLst.txt` specifies the mappable, non-blacklisted genome we use for randomization analyses

### Code

#### Structural Variant Merging
- `script_merge_translocations[...].py` collapses interchromosomal structural variants from different callers
- `script_filter_mapqBAM_3.0.py` filters collapsed interchromosomal structural variants using common criteria
- `script_merge_intraSVs_[...].py` collapses intrachromosomal structural variants from different callers
- `script_filter_atleast2callers_2.0.py` filters collapsed intrachromosomal structural variants using common criteria
- `script_merge_circularized_regions_[...].py` overlaps interchromosomal structural variants with circle calls and classifies accordingly
- `script_merge_circularized_regions_[...]_intraSV_[...].py` overlaps intrachromosomal structural variants with circle calls and classifies accordingly
- `script_search_integrations_[...].py` searches for two-breakpoint patterns indicative of circle integration into a chromosome

#### Detection of Clusters of Rearrangements
Tree-shaped rearrangement regions are referred to as *Palm Trees* throughout the code.
- `RunAll.R` runs scripts in the correct order
- `Parse*.R` and `Prepare*.R` are several scripts to read data from different sources and create tidy representations for further analysis
- `CallPalmTrees.R` uses output of structural variant merging to call palm trees (a.k.a. tree-shaped rearrangements or clusters of rearrangements)
- `CallPalmTreesPedPanCan.R` calls palm trees on published data from Gröbner et al. 2018 (using the same settings as in CallPalmTrees.R)


#### Basic Rearrangement Cluster Statistics
- `CircosPlots.R` plots the identified rearrangements and palm tree regions
- `PalmTreeStackPlot.R` creates plots that integrate all CN, SV and Circle-seq information
- `GeneralPalmTreeStatistics.R` analyses number and length distributions of palm trees
- `PalmTreeDensity.R` plots genome-wide recurrence of palm trees
- `PedPanCanStats.R` explores palm tree occurrence in PedPanCan dataset
- `CompareBerlinCohortToPedPanCan.R` compares palm tree prevalence between the Peifer/Berlin dataset and the PedPanCan neuroblastoma cases
- `AnalyseNBCallPalmTreesRandomization.Rmd` explores 500 synthetic datasets for our cohort and the PedPanCan dataset respectively. Synthetic datasets were obtained by randomizing breakpoint positions during palm tree calling. This is used to estimate false-discovery rates due to the number of rearrangements per sample.

#### Genes and Gene Expression
- `PalmTreeGenes.R` analyzes which genes can be found within palm tree regions
- `PalmTreeTargetsCloseToOncogenes.R` statistically analyses whether palm tree target sites are significantly enriched in the neighborhood of cancer-related genes.
- `ExpressionAnalysis.R` screens for deregulated genes close to palm tree associated rearrangements

#### Association with circular DNA
- `RegionSampling.R` contains sampling methods to randomly sample regions from a masked or unmasked genome
- `Overlap[...]CircleCalls.R` statistically analyses the degree of overlap between palm tree regions and eccDNA / ecDNA inferred from WGS / Circle-seq
- `MergeWGSCSOverlapPlots.R` merges results from the preceding scripts for integrative plotting

#### Clinical Relevance
- `ClinicalData.R` explores the prognostic significance of palm trees

#### Breakpoint Analysis
- `CircleJunctionAnalyseSvabaBreakpoints.R` investigates sequence characteristics reconstructed Circle-seq circle junctions
- `SVAnalyseSvabaBreakpoints.R` investigates sequence characteristics interchromosomal rearrangements calls by Svaba
- `MemeAnalysis.sh` runs MEME on several breakpoint-associated sequences
- `myHomology.R` obtains microhomology estimates for accurate breakpoint coordinates from the reference genome
- `Repeats.R` tests for association of breakpoints with repetitive regions

### Contact
If you have any questions concerning code or data, please do not hesitate to contact us at henssenlab@gmail.com.
