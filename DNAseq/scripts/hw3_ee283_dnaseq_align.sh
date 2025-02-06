#!/bin/bash
#SBATCH --job-name=align_DNAseq    ## Name of the job
#SBATCH -A CLASS-ECOEVO283         ## Account to charge
#SBATCH -p standard                ## Partition/queue name
#SBATCH --array=1-12               ## Number of tasks to launch (matches number of lines in prefixes_DNAseq.txt)
#SBATCH --cpus-per-task=2          ## Number of cores the job needs
#SBATCH --error=logs/dna_alignment.%A_%a.err
#SBATCH --output=logs/dna_alignment.%A_%a.out
#SBATCH --time=02:00:00            ## Time limit
#SBATCH --mem=4G                   ## Memory required per node

module load bwa/0.7.8
module load samtools/1.10

# Define paths
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
rawDataDir="/pub/$USER/EE283/DNAseq/rawdata"
outputDir="/pub/$USER/EE283/DNAseq/output"

# Ensure output directory exists
mkdir -p $outputDir

# Read the prefix for the current array task
file="/pub/$USER/EE283/DNAseq/scripts/prefixes_DNAseq.txt"
prefix=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file)

# Align data with BWA MEM, using multiple threads and marking secondary alignments
bwa mem -t $SLURM_CPUS_PER_TASK -M $ref $rawDataDir/${prefix}_1.fq.gz $rawDataDir/${prefix}_2.fq.gz | \
samtools view -bS - > $outputDir/${prefix}.bam

# Sort the BAM file
samtools sort $outputDir/${prefix}.bam -o $outputDir/${prefix}.sorted.bam

# Index the sorted BAM file
samtools index $outputDir/${prefix}.sorted.bam

# Clean up intermediate file
rm $outputDir/${prefix}.bam