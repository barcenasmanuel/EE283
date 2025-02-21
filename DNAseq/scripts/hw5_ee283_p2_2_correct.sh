#!/bin/bash
#SBATCH --job-name=merge_advanced_filtered_vcfs
#SBATCH --account=CLASS-ECOEVO283
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/merge_advanced_filtered_vcfs_sam_%A.out
#SBATCH --error=logs/merge_advanced_filtered_vcfs_sam_%A.err

module load bcftools/1.15.1
module load samtools/1.15.1

outputDir="/pub/$USER/EE283/DNAseq/output/snp_results"
advancedFilteredDir="${outputDir}/advanced_filtered_snps"
chromNamesFile="${outputDir}/chrome.names.txt"

# Merge all advanced filtered VCFs
bcftools concat -O v $(for chr in $(cat $chromNamesFile); do echo -n "${advancedFilteredDir}/${chr}_final_filtered.vcf.gz "; done) | \
bgzip -c > ${advancedFilteredDir}/all_chromosomes_advanced_filtered.vcf.gz

# Index the merged VCF file
tabix -p vcf ${advancedFilteredDir}/all_chromosomes_advanced_filtered.vcf.gz

echo "Final merged advanced filtered VCF: ${advancedFilteredDir}/all_chromosomes_advanced_filtered.vcf.gz"
echo "Index file created: ${advancedFilteredDir}/all_chromosomes_advanced_filtered.vcf.gz.tbi"

# Count SNPs in the merged file
final_count=$(bcftools view -H ${advancedFilteredDir}/all_chromosomes_advanced_filtered.vcf.gz | wc -l)
echo "Total SNP count in merged file: $final_count"