#!/bin/bash
#SBATCH --job-name=index_genome_rna
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/index_genome.out
#SBATCH --error=logs/index_genome.err

module load hisat2/2.2.1
module load python/3.10.2

# Define paths to your reference genome and GTF file
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
gtf="/pub/$USER/EE283/ref/dmel-all-r6.13.gtf"

# Extract splice sites and exons
python hisat2_extract_splice_sites.py $gtf > ref/dmel.ss
python hisat2_extract_exons.py $gtf > ref/dmel.exon

# Build the Hisat2 index
hisat2-build -p 8 --exon ref/dmel.exon --ss ref/dmel.ss $ref ref/dmel_trans