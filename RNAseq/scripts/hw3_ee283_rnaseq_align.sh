#!/bin/bash
#SBATCH --job-name=align_RNAseq
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --array=1-20                ## Number of tasks to launch (10 samples x 2 tissues)
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/rnaseq_alignment.%A_%a.out
#SBATCH --error=logs/rnaseq_alignment.%A_%a.err

module load hisat2/2.2.1
module load samtools/1.10

# Define paths
index="ref/dmel_trans"
rawDataDir="/pub/$USER/EE283/RNAseq/rawdata"
outputDir="/pub/$USER/EE283/RNAseq/output"

# Ensure output directory exists
mkdir -p $outputDir

# Sample IDs to process
samples=("21148" "21286" "22162" "21297" "21029" "22052" "22031" "21293" "22378" "22390")

# Read the sample ID and tissue type for the current array task
sample_index=$((($SLURM_ARRAY_TASK_ID - 1) / 2))
tissue_index=$((($SLURM_ARRAY_TASK_ID - 1) % 2))
sample_id=${samples[$sample_index]}
tissue_type=("E" "B")

# Define input files based on tissue type
input_R1="${rawDataDir}/${sample_id}_${tissue_type[$tissue_index]}_0_R1.fastq.gz"
input_R2="${rawDataDir}/${sample_id}_${tissue_type[$tissue_index]}_0_R2.fastq.gz"

# Align reads with Hisat2
hisat2 -p $SLURM_CPUS_PER_TASK -x $index -1 $input_R1 -2 $input_R2 | \
samtools view -bS - > $outputDir/${sample_id}_${tissue_type[$tissue_index]}.bam

# Sort the BAM file
samtools sort $outputDir/${sample_id}_${tissue_type[$tissue_index]}.bam -o $outputDir/${sample_id}_${tissue_type[$tissue_index]}.sorted.bam

# Index the sorted BAM file
samtools index $outputDir/${sample_id}_${tissue_type[$tissue_index]}.sorted.bam

# Clean up intermediate file
rm $outputDir/${sample_id}_${tissue_type[$tissue_index]}.bam