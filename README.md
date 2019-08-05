This repository contains code and data to reproduce the findings reported on tree-shaped rearrangements in Koche et al. (under review).

### Data

*Data from the Berlin Neuroblastoma Dataset can only be shared upon request at the moment. This is due to ongoing work in collaborating groups. It will be made public at the time of publication.*

**Circle Calls**
- `Data/CircleseqCircles.bed` contains Circle-seq circle calls for our dataset
- `Data/WGSinferredCircles_eccDNA.txt` contains WGS-inferred circle calls for our dataset (eccDNA only)
- `Data/WGSinferredCircles_ecDNA.txt` contains WGS-inferred circle calls for our dataset (ecDNA only)
- `Data/CircleSeq_CircleReads_Svaba/` contains Svaba-based SV calls based on circle-supporting reads for base pair-accurate reassembly of circle junctions

**WGS Rearrangement Calls**
- `Data/MergedSV_StrictestFiltering_CircleSeqCircleAnnotation.txt` are the merged and filtered interchromosomal rearrangement calls (see Structural Variant Merging) which were overlapped with Circle-seq circle calls for classification *(available upon request)*
- `Data/MergedSV_StrictestFiltering_WGSCircleAnnotation.txt` are the merged and filtered interchromosomal rearrangement calls (see Structural Variant Merging) which were overlapped with WGS-inferred circle calls for classification *(available upon request)*
- `Data/PedPanCanSVs.csv` contains structural variants in the DKFZ Pediatric Pan-Cancer dataset published by Gröbner et al. 2018 (data downloaded from the corresponding [R2 platform](https://hgserver1.amc.nl/cgi-bin/r2/main.cgi?&dscope=DKFZ_PED&option=about_dscope))
- `Data/PedPanCanMeta.csv` contains metadata such as the tumor entity for the DKFZ Pediatric Pan-Cancer dataset

**Expression Data**
- `Data/berlin_cohort_rnaseq_fpkm.txt` is RNA-seq data for the Berlin neuroblastoma cohort *(available upon request)*
- `Data/peifer_54nb_fpkms.txt` is RNA-seq data published by Peifer et al.

**Clinical Data**
- `Data/ClinicalData.csv` specifies risk group and survival data for each patient

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
- `RunAll.R` runs all required scripts in the correct order
- `Parse*.R` are several scripts to read data from different sources and create tidy representations thereof
- `CallPalmTrees.R` uses output of structural variant merging to call palm trees (a.k.a. clusters of rearrangements)
- `CallPalmTreesPedPanCan.R` calls palm trees on published data from Gröbner et al. 2018 (using the same settings as in CallPalmTrees.R)
- `AnalyseNBCallPalmTreesRandomization.Rmd` explores 500 synthetic datasets for our cohort and the PedPanCan dataset respectively. Synthetic datasets were obtained by randomizing breakpoint positions during palm tree calling. This is used to estimate false-discovery rates due to the number of rearrangements per sample.

#### Basic Rearrangement Cluster Statistics
- `CircosPlots.R` plots the identified rearrangements and palm tree regions
- `PalmTreeStackPlot.R` creates plots that integrate all CN, SV and Circle-seq information
- `GeneralPalmTreeStatistics.R` analyses number and length distributions of palm trees
- `PalmTreeDensity.R` plots genome-wide recurrence of palm trees

#### Genes and Gene Expression
- `PalmTreeGenes.R` analyzes which genes can be found within palm tree regions
- `PalmTreeTargetsCloseToOncogenes.R` statistically analyses whether palm tree target sites are significantly enriched in the neighborhood of cancer-related genes.
- `ExpressionAnalysis.R` screens for deregulated genes close to palm tree associated rearrangements

#### Association with circular DNA
- `RegionSampling.R` contains sampling methods to randomly sample regions from masked or unmasked genome
- `OverlapCircleSeqCircleCalls.R` statistically analyses the degree of overlap between palm tree regions and circles as called from Circle-Seq
- `OverlapWGSCircleCalls.R` statistically analyses the degree of overlap between palm tree regions and circles as called from WGS
- `MergeWGSCSOverlapPlots.R` merges results from the two preceding scripts for integrative plotting

#### Clinical Relevance
- `ClinicalData.R` explores the prognostic significance of palm trees

#### Breakpoint Analysis
- `SVAnalyseSvabaBreakpoints.R` investigates sequence characteristics of base pair-accurate reconstruction of Circle-seq circle junctions
- `MemeAnalysis.sh` runs MEME on several breakpoint-associated sequences
- `myHomology.R` obtains microhomology estimates for accurate breakpoint coordinates from the reference genome

### Contact
If you have any questions concerning code or data, please do not hesitate to contact us at henssenlab@gmail.com.
