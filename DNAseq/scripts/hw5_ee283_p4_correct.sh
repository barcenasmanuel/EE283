#!/bin/bash
#SBATCH --job-name=extract_filter_snps_task4_advanced
#SBATCH --account=CLASS-ECOEVO283
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/extract_filter_snps_task4_advanced_%A.out
#SBATCH --error=logs/extract_filter_snps_task4_advanced_%A.err

module load vcftools/0.1.16
module load bcftools/1.15.1

# Define paths
inputDir="/pub/$USER/EE283/DNAseq/output/snp_results/advanced_filtered_snps"
outputDir="/pub/$USER/EE283/DNAseq/output/snp_results/visualization_task4_advanced"
mkdir -p $outputDir

# Input VCF file (X chromosome, advanced filtered)
input_vcf="${inputDir}/X_final_filtered.vcf.gz"

# Extract biallelic SNPs from the first 1Mb of the X chromosome
bcftools view -r X:1-1000000 -m2 -M2 -v snps -O z -o ${outputDir}/X_1Mb_task4.vcf.gz $input_vcf

# Convert to 012 format
vcftools --gzvcf ${outputDir}/X_1Mb_task4.vcf.gz \
    --012 \
    --out ${outputDir}/X_1Mb_012_task4

# Count SNPs before filtering
total_snps=$(bcftools view -H ${outputDir}/X_1Mb_task4.vcf.gz | wc -l)
echo "Total SNPs in the first 1Mb of X chromosome: $total_snps"

# Print the first few lines of the 012 file
echo "First few lines of the 012 file:"
head -n 5 ${outputDir}/X_1Mb_012_task4.012

# Count SNPs with no heterozygotes
no_het_count=$(awk '($3!=1 && $4!=1 && $5!=1 && $6!=1)' ${outputDir}/X_1Mb_012_task4.012 | wc -l)
echo "SNPs with no heterozygotes: $no_het_count"

# Count SNPs with sum of genotypes equal to 4
sum_4_count=$(awk '($3+$4+$5+$6==4)' ${outputDir}/X_1Mb_012_task4.012 | wc -l)
echo "SNPs with sum of genotypes equal to 4: $sum_4_count"

# Filter SNPs using modified awk command
awk 'NR==1 || (($3==0 || $3==2) && ($4==0 || $4==2) && ($5==0 || $5==2) && ($6==0 || $6==2)) && 
    (($3+$4+$5+$6==4) || ($3+$4+$5+$6==2 && ($3==-1 || $4==-1 || $5==-1 || $6==-1)))' \
    ${outputDir}/X_1Mb_012_task4.012 > ${outputDir}/X_1Mb_012_filtered_task4.012

# Count SNPs after new filtering
new_filtered_snps=$(wc -l < ${outputDir}/X_1Mb_012_filtered_task4.012)
new_filtered_snps=$((new_filtered_snps - 1))  # Subtract 1 to account for the header line

echo "SNPs after new filtering criteria: $new_filtered_snps"

# Print the filtered SNPs
if [ $new_filtered_snps -gt 0 ]; then
    echo "Filtered SNPs:"
    cat ${outputDir}/X_1Mb_012_filtered_task4.012
else
    echo "No SNPs met the new filtering criteria."
fi

echo "Task 4 Advanced: Data extraction, formatting, and filtering complete. Output files are in $outputDir"