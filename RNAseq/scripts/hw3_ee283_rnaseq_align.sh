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
index="/pub/$USER/EE283/ref/dmel_trans"
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

# Find the appropriate replicate files
input_R1=$(find $rawDataDir -name "${sample_id}_${tissue_type[$tissue_index]}_*_R1.fastq.gz" | head -n 1)
input_R2=$(find $rawDataDir -name "${sample_id}_${tissue_type[$tissue_index]}_*_R2.fastq.gz" | head -n 1)

# Check if files were found
if [[ -z $input_R1 || -z $input_R2 ]]; then
  echo "Error: Input files for ${sample_id}_${tissue_type[$tissue_index]} not found."
  exit 1
fi

# Extract replicate number from the filename
replicate=$(basename $input_R1 | cut -d'_' -f3)

# Define output file names
bam_file="${outputDir}/${sample_id}_${tissue_type[$tissue_index]}_${replicate}.bam"
sorted_bam_file="${outputDir}/${sample_id}_${tissue_type[$tissue_index]}_${replicate}.sorted.bam"

# Remove old output files if they exist
rm -f $bam_file $sorted_bam_file $sorted_bam_file.bai

# Align reads with Hisat2
hisat2 -p $SLURM_CPUS_PER_TASK -x $index -1 $input_R1 -2 $input_R2 | \
samtools view -bS - > $bam_file

# Sort the BAM file
samtools sort $bam_file -o $sorted_bam_file

# Index the sorted BAM file
samtools index $sorted_bam_file

# Clean up intermediate file
rm $bam_file