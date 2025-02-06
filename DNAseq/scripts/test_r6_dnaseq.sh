#!/bin/bash
#SBATCH --job-name=test_alignment    ## Name of the job
#SBATCH -A CLASS-ECOEVO283           ## Account to charge
#SBATCH -p standard                  ## Partition/queue name
#SBATCH --cpus-per-task=2            ## Number of cores the job needs
#SBATCH --error=test_alignment.%J.err
#SBATCH --output=test_alignment.%J.out
#SBATCH --time=01:00:00              ## Time limit
#SBATCH --mem=4G                     ## Memory required per node

module load bwa/0.7.8
module load samtools/1.10

# Define paths
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
rawDataDir="/pub/$USER/EE283/DNAseq/rawdata"
outputDir="/pub/$USER/EE283/DNAseq/test_output"

# Ensure output directory exists
mkdir -p $outputDir

# Align test data with BWA MEM, using multiple threads and marking secondary alignments
bwa mem -t $SLURM_CPUS_PER_TASK -M $ref $outputDir/test_1.fq.gz $outputDir/test_2.fq.gz | \
samtools view -bS - > $outputDir/test.bam

# Sort the BAM file
samtools sort $outputDir/test.bam -o $outputDir/test.sorted.bam

# Clean up intermediate file
rm $outputDir/test.bam

# Optional: sleep for 2 minutes
sleep 2m