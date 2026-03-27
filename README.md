# ChloroplastRedoxReg

RNA-Seq analysis of Arabidopsis leaves illuminated with visible monochromatic light and combinations of visible wavelengths forming white light favoring either Photosystem II (PSII) or Photosystem I (PSI), where we studied the transcriptional regulation by the redox state of the plastoquinone (PQ) pool. 

The repository includes both standalone **R scripts** and a **modular Nextflow pipeline** for scalable, reproducible computation on HPC systems or local environments using Docker or Singularity.

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

---------------

# Nextflow RNA‑Seq analysis workflow

A fully modular and reproducible workflow implemented in **Nextflow DSL2**.

## Pipeline modules 

### **1. LOAD_FILTER**
- Loads `.RData` files  
- Filters low-expression genes  
- Outputs `dds1h.rds` and `dds4w.rds`

### **2. DESEQ_PAIRWISE**
- Performs pairwise DESeq2 contrasts  
- LFC shrinkage with apeglm  
- Outputs CSV files of significant genes

### **3. DESEQ_GROUP**
- Performs group-level comparisons (g1 vs g2)
- Outputs grouped DE genes

### **4. GO_ENRICHMENT**
- Performs GO enrichment for up- and downregulated genes  
- Outputs GO term tables

---

# Running the pipeline

### With Docker

nextflow run main.nf -with-docker deseq2-pipeline

### Without containers (local R installation)

nextflow run main.nf

---

# Docker environment

The provided Dockerfile uses **rocker/tidyverse** and installs:

- DESeq2  
- apeglm  
- GSEABase  
- Category  
- GOstats

Build locally with:

docker build -t deseq2-pipeline .
