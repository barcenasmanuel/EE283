#!/bin/bash
#SBATCH --job-name=align_ATACseq    ## Name of the job
#SBATCH -A CLASS-ECOEVO283          ## Account to charge
#SBATCH -p standard                 ## Partition/queue name
#SBATCH --array=1-24                ## Number of tasks to launch (matches number of lines in prefixes_ATACseq.txt)
#SBATCH --cpus-per-task=2           ## Number of cores the job needs
#SBATCH --error=logs/atacseq_alignment.%A_%a.err
#SBATCH --output=logs/atacseq_alignment.%A_%a.out
#SBATCH --time=02:00:00             ## Time limit
#SBATCH --mem=4G                    ## Memory required per node

module load bwa/0.7.8
module load samtools/1.10

# Define paths
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
rawDataDir="/pub/$USER/EE283/ATACseq/rawdata"
outputDir="/pub/$USER/EE283/ATACseq/output"

# Ensure output directory exists
mkdir -p $outputDir

# Read the prefix for the current array task
file="/pub/$USER/EE283/ATACseq/scripts/prefixes_ATACseq.txt"
prefix=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $file)

# Align reads directly with BWA MEM
bwa mem -t $SLURM_CPUS_PER_TASK -M $ref $rawDataDir/${prefix}_R1.fq.gz $rawDataDir/${prefix}_R2.fq.gz | \
samtools view -bS - > $outputDir/${prefix}.bam

# Sort the BAM file
samtools sort $outputDir/${prefix}.bam -o $outputDir/${prefix}.sorted.bam

# Index the sorted BAM file
samtools index $outputDir/${prefix}.sorted.bam

# Clean up intermediate file
rm $outputDir/${prefix}.bam