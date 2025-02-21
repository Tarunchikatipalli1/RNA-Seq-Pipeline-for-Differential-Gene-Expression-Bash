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
### Steps in `RNAseq_pipeline.sh:`
- FastQC: Checks raw .fastq quality, outputs to output `/1_initial_qc/`.
- Trimmomatic: Trims reads with a Phred score < 10, outputs to `output/2_trimmed_output/`.
- HISAT2: Aligns trimmed reads to GRCm38, converts SAM to sorted BAM, outputs to `output/3_aligned_sequences/aligned_bam/`.
- featureCounts: Quantifies gene counts, outputs to `output/4_final_counts/final_counts.txt`.
- MultiQC: Generates a report in `output/5_multiQC/`.
  
### Notes:
- Update paths in the script (e.g., `/Users/tarun/Desktop/RNASeq_pipeline/`) to match your system.
- Ensure all tools are in your PATH or provide full paths.

### 2. Differential Gene Expression Analysis
Run the R script in RStudio to analyze gene counts:

```bash
source("RNAseq_analysis.R")
```
### Steps in `RNAseq_analysis.R`:
- **Setup**: Installs and loads required R packages.
- **Data Import**: Loads `final_counts.txt` and `metadata.txt` (update paths as needed).
- **DESeq2**: Performs DGE analysis comparing groups (e.g., LoGlu vs. HiGlu).
- **Annotation**: Adds gene names, Entrez IDs, and Ensembl IDs using `org.Mm.eg.db`.
- **Output**: Writes normalized counts and annotated results to `.txt` files.
- **Visualization**: Generates PCA, heatmap, volcano, MA, dispersion, and single-gene plots.
- **Pathway Analysis**: Performs KEGG and GO enrichment on significant genes.
  
**Notes:**
- Replace `/Users/tarun/Desktop/RNASeq_pipeline/example/` with your data paths.
- Plots are displayed in RStudio; save them manually if needed (e.g., `ggsave()`).

## Output
### Bash Pipeline:
- `output/1_initial_qc/`: FastQC HTML and ZIP files.
- `output/2_trimmed_output/`: Trimmed `.fq` files.
- `output/3_aligned_sequences/aligned_bam/`: Sorted `.bam` files.
- `output/4_final_counts/final_counts.txt`: Gene count table.
- `output/5_multiQC/multiqc_report.html`: QC report.

### R Analysis:
- `normalized_counts.txt`: Normalized gene counts.
- `normalized_counts_significant.txt`: Significant gene counts.
- `results_gene_annotated.txt`: Full annotated results.
- `results_gene_annotated_significant.txt`: Significant annotated results.
- Plots: PCA, heatmap, volcano, etc. (view in RStudio).

## Customization
- **Input Data**: Place your `.fastq` files in `input/` or update the script.
- **Genome**: Use a different HISAT2 index for other organisms (e.g., human GRCh38).
- **Parameters**: Adjust Trimmomatic (`TRAILING:10`) or DESeq2 (`alpha = 0.05`) settings as needed.










