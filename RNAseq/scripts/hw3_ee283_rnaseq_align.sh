#!/bin/bash
#SBATCH --job-name=align_RNAseq
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --array=1-20                ## Number of tasks to launch (for the 20 samples)
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/rnaseq_alignment.%A_%a.out
#SBATCH --error=logs/rnaseq_alignment.%A_%a.err

module load hisat2/2.2.1
module load samtools/1.10

# Define paths
index="ref/dm6_trans"
rawDataDir="/pub/$USER/EE283/RNAseq/rawdata"
outputDir="/pub/$USER/EE283/RNAseq/output"

# Ensure output directory exists
mkdir -p $outputDir

# Sample IDs to process
samples=("21148" "21286" "22162" "21297" "21029" "22052" "22031" "21293" "22378" "22390")

# Read the sample ID for the current array task
sample_id=${samples[$SLURM_ARRAY_TASK_ID-1]}

# Define input files
input_R1="${rawDataDir}/${sample_id}_E_R1.fq.gz"
input_R2="${rawDataDir}/${sample_id}_E_R2.fq.gz"

# Align reads with Hisat2
hisat2 -p $SLURM_CPUS_PER_TASK -x $index -1 $input_R1 -2 $input_R2 | \
samtools view -bS - > $outputDir/${sample_id}.bam

# Sort the BAM file
samtools sort $outputDir/${sample_id}.bam -o $outputDir/${sample_id}.sorted.bam

# Index the sorted BAM file
samtools index $outputDir/${sample_id}.sorted.bam

# Clean up intermediate file
rm $outputDir/${sample_id}.bam