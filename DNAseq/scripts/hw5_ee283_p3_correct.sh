#!/bin/bash
#SBATCH --job-name=extract_snps_advanced
#SBATCH --account=CLASS-ECOEVO283
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/extract_snps_advanced_%A.out
#SBATCH --error=logs/extract_snps_advanced_%A.err

module load vcftools/0.1.16
module load bcftools/1.15.1

# Define paths
inputDir="/pub/$USER/EE283/DNAseq/output/snp_results/advanced_filtered_snps"
outputDir="/pub/$USER/EE283/DNAseq/output/snp_results/visualization_advanced"
mkdir -p $outputDir

# Input VCF file (X chromosome, advanced filtered)
input_vcf="${inputDir}/X_final_filtered.vcf.gz"

# Extract biallelic SNPs from the first 1Mb of the X chromosome
bcftools view -r X:1-1000000 -m2 -M2 -v snps -O z -o ${outputDir}/X_1Mb.vcf.gz $input_vcf

# Convert to 012 format
vcftools --gzvcf ${outputDir}/X_1Mb.vcf.gz \
    --012 \
    --out ${outputDir}/X_1Mb_012

echo "Data extraction and formatting complete. Output files are in $outputDir"

# Count SNPs
snp_count=$(bcftools view -H ${outputDir}/X_1Mb.vcf.gz | wc -l)
echo "Number of SNPs in the first 1Mb of X chromosome: $snp_count"