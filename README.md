# ChloroplastRedoxReg

# RNA-Seq analysis of Arabidopsis leaves illuminated with visible monochromatic light and combinations of visible wavelengths forming white light favoring either Photosystem II (PSII) or Photosystem I (PSI), where we studied the transcriptional regulation by the redox state of the plastoquinone (PQ) pool. 

This repository contains a collection of scripts for RNA-seq data analysis, including transcript quantification, differential expression, and data visualization. The analysis pipeline makes use of high-performance computing (HPC) clusters and R scripts.

## Scripts Overview

### 1. `STAR_quantTranscript.sh` Script
- **Description**: This script utilizes the HPC cluster at the Center for Scientific Computing (Finland) to process and quantify sequencing files in `fastq` format.
- **Workflow**:
  - Generates an index using the **STAR** aligner.
  - Aligns reads to the genome index to produce SAM output.
  - Uses **RSEM** to quantify transcripts based on the aligned reads.
  
### 2. `cluster_profiles.R`
- **Description**: This R script plots the expression of genes for each identified cluster separately. It helps visualize gene expression across different gene clusters.

### 3. `difexpr.R`
- **Description**: Performs differential expression analysis using **DESeq2**. This script identifies differentially expressed genes between sample groups and prepares the results for downstream analysis.

### 4. `Figures/`
- **Description**: Contains figures generated from the differential expression and splicing analysis. These include visualizations of significant genes, cluster profiles, and other relevant results.
