#!/bin/bash

# Enable exit on error
set -e

# Start timing
SECONDS=0


# Step 1: Quality check with FastQC
echo "Running FastQC..."
fastqc -o /Users/tarun/Desktop/RNASeq_pipeline/output/1_initial_qc/ --noextract /Users/tarun/Desktop/RNASeq_pipeline/input/*.fastq
echo "FastQC completed."

# Step 2: Trim low-quality sequences with Trimmomatic
echo "Running Trimmomatic..."
for file in input/*.fastq; do
    base=$(basename "$file" .fastq)
    java -jar /Users/tarun/Desktop/Tools/trimmomatic-0.39.jar SE -threads 4 \
        "$file" /Users/tarun/Desktop/RNASeq_pipeline/output/2_trimmed_output/"${base}_trimmed.fq" \
        TRAILING:10 -phred33
    echo "Trimming completed for $file."
done
echo "Trimmomatic processing complete!"

# Step 3: Run HISAT2 alignment
echo "Running HISAT2 alignment with pre-built index..."
for file in output/2_trimmed_output/*.fq; do
    base=$(basename "$file" .fq)
    
    # Run HISAT2 alignment to generate SAM files
    hisat2 -p 4 -x /Users/tarun/Desktop/RNASeq_pipeline/HISAT2/grcm38/genome -U "$file" -S /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}.sam"
    echo "HISAT2 alignment completed for $file."

    # Convert SAM to BAM
    samtools view -Sb /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}.sam" > /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}.bam"
    echo "SAM to BAM conversion completed for ${base}.sam."

    # Sort the BAM file
    samtools sort /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}.bam" -o /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}_sorted.bam"
    echo "BAM sorting completed for ${base}.bam."

    # Remove intermediate SAM and unsorted BAM files to save space
    rm /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}.sam" /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/"${base}.bam"
    echo "Removed intermediate files for ${base}."

done
echo "HISAT2 alignment and BAM processing completed."

# Step 4: Gene counting with featureCounts
echo "Running featureCounts..."
featureCounts -a /Users/tarun/Desktop/RNASeq_pipeline/annotation/gencode.vM12.annotation.gtf -o /Users/tarun/Desktop/RNASeq_pipeline/output/4_final_counts/final_counts.txt -T 4 /Users/tarun/Desktop/RNASeq_pipeline/output/3_aligned_sequences/aligned_bam/*_sorted.bam
echo "featureCounts completed."

# Step 5: Generate MultiQC report
echo "Generating MultiQC report..."
multiqc output --outdir /Users/tarun/Desktop/RNASeq_pipeline/output/5_multiQC
echo "MultiQC report generated."

# Display elapsed time
duration=$SECONDS
echo "RNA-seq pipeline completed successfully!"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."