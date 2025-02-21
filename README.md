# RNA-seq Pipeline and Analysis


This repository contains my custom RNA-seq pipeline and analysis workflow, designed to process RNA-seq data and perform differential gene expression (DGE) analysis. It includes a Bash script for preprocessing and alignment, and an R script for statistical analysis and visualization.

## Overview
This project processes RNA-seq data through the following steps:  
1. **Quality Control**: FastQC for raw sequence quality checking.  
2. **Trimming**: Trimmomatic to remove low-quality reads.  
3. **Alignment**: HISAT2 to align reads to the GRCm38 mouse genome.  
4. **Quantification**: featureCounts for gene counting.  
5. **Reporting**: MultiQC for summary reports.  
6. **DGE Analysis**: DESeq2 in R for differential expression, followed by visualization and pathway analysis.  

The pipeline is configured for my local setup (`/Users/tarun/Desktop/RNASeq_pipeline/`) but can be adapted by modifying file paths.

## Repository Contents
- **`RNAseq_pipeline.sh`**: Bash script for preprocessing, alignment, and quantification.  
- **`RNAseq_analysis.R`**: R script for DGE analysis, visualization, and pathway enrichment.  
- **`example/`**: Directory with sample `final_counts.txt` and `metadata.txt` files for testing the R script (replace with your actual data paths).  

## Prerequisites

### Software
- **Bash Tools**:  
  - `fastqc`  
  - `trimmomatic-0.39.jar`  
  - `hisat2`  
  - `samtools`  
  - `featureCounts`  
  - `multiqc`
    
  (Install via Bioconda: `conda install -c bioconda fastqc trimmomatic hisat2 samtools subread multiqc`)  
- **R**: Install R and RStudio (https://www.rstudio.com/).  
- **R Packages**: Listed in `RNAseq_analysis.R` (e.g., DESeq2, ggplot2, clusterProfiler, etc.). Install them by running the script’s setup section.

### Data
- **Input**: `.fastq` files in `input/` (e.g., from SRA or your sequencing run).  
- **Reference**: GRCm38 genome index for HISAT2 (`HISAT2/grcm38/genome`).  
- **Annotation**: Gencode vM12 GTF (`annotation/gencode.vM12.annotation.gtf`).  

## Usage

### 1. Preprocessing and Alignment
Run the Bash script to process your RNA-seq data:  

```bash
bash RNAseq_pipeline.sh
```
