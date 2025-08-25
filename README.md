# ISR_GDF15_project

Analysis of Integrated Stress Response (ISR) and GDF15 gene expression patterns across fibroblast lifespan and GTEx tissue data.

## Overview

This repository contains R scripts and analysis for investigating:
- **Fibroblast lifespan**: RNA-seq analysis of cellular lifespan and ISR gene expression
- **GTEx tissue comparison**: Age-related gene expression patterns across human tissues
- **ISR pathway analysis**: Factor analysis and correlation studies of stress response genes

## Structure

- `Data/` - Input datasets include the Fibroblast lifespan RNA-seq, Opto-PKR RNA-seq, and GTEx tissue data RNA-seq. Also contains Processed data from some analyses required for later analyses. 
- `Scripts/` - R scripts for data processing and analysis
- `Results/` - Output files, plots, and analysis results. 

## Key Analyses

- ISR-GDF15 index generation: Factor analysis of ISR gene expression cellular lifespan studies with various stress conditions
- Opto-PKR validation: Validation of the fibroblast-derived ISR index using optogenetic PKR activation data
- ISR-GDF15 index applied across human tissues
- Spearman correlation analysis between fibroblast and GTEx aging patterns
- Age-related gene expression comparisons across human tissues


## Requirements

- R with tidyverse, edgeR, ComplexHeatmap, and other bioconductor packages
- See individual scripts for specific package dependencies

REQUIRED DOWNLOADS!!!
- The fibroblast lifespan data needs to be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179848 and put into the folder named Fibroblast_lifespan in Data. The file to download is: GSE179848_processed_cell_lifespan_RNAseq_data.csv.gz
- The GTEx RNA-seq data can be acquired from the GTEx Portal https://www.gtexportal.org/home/ Individual tissues were put into their own tissue folders for the code to run smoothly. This must be put into the All_Tissues_Indiv_Folders folder in the gtex folder within the Data folder.
- Access to Protected Access Data must be applied for at https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424.v10.p2 This must be put into the gtex folder in the Data folder.
             
	

## Usage

Run scripts from the `Scripts/` directory. Main analysis scripts:
- `Fibroblast_lifespan/` - Cellular lifespan and ISR expression studies
- `Opto_PKR/` - Validation data using optogenetic PKR activation
- `gtex/` - Exploration of ISR-GDF15 index used on human tissue RNA-seq, investigating differences in ISR-GDF15 score based on place of death, and proliferation status of tisssues
- `AgingRhos_for_correlation_FB_GTEX.R` - GTEx vs fibroblast correlation analysis of ISR-GDF15 correlations with age


