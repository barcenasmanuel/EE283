#!/bin/bash
#SBATCH --job-name=advanced_filter_snps
#SBATCH --account=CLASS-ECOEVO283
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=1-7
#SBATCH --output=logs/advanced_filter_snps_index_sam_%A_%a.out
#SBATCH --error=logs/advanced_filter_snps_index_sam_%A_%a.err

module load bcftools/1.15.1
module load samtools/1.15.1

# Define paths
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
inputDir="/pub/$USER/EE283/DNAseq/output"
outputDir="/pub/$USER/EE283/DNAseq/output/snp_results"
chromNamesFile="${outputDir}/chrome.names.txt"

# Create a directory for advanced filtered SNPs
advancedFilteredDir="${outputDir}/advanced_filtered_snps"
mkdir -p $advancedFilteredDir

# Get the chromosome for this array task
chr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chromNamesFile)
input_vcf="${outputDir}/result.${chr}.vcf.gz"

# Count initial SNPs
initial_count=$(bcftools view -H $input_vcf | wc -l)
echo "Initial SNP count for ${chr}: $initial_count"

# Step 1: Apply advanced filters
bcftools filter -i 'FS<40.0 && SOR<3 && MQ>40.0 && MQRankSum>-5.0 && MQRankSum<5 && QD>2.0 && ReadPosRankSum>-4.0 && INFO/DP<16000' \
    -O v ${input_vcf} | bgzip -c > ${advancedFilteredDir}/${chr}_step1_advanced_filtered.vcf.gz
tabix -p vcf ${advancedFilteredDir}/${chr}_step1_advanced_filtered.vcf.gz

step1_count=$(bcftools view -H ${advancedFilteredDir}/${chr}_step1_advanced_filtered.vcf.gz | wc -l)
echo "SNP count after advanced filtering for ${chr}: $step1_count"

# Step 2: Filter SNPs with low depth or genotype quality in some individuals
bcftools filter -S . -e 'FMT/DP<3 | FMT/GQ<20' \
    -O v ${advancedFilteredDir}/${chr}_step1_advanced_filtered.vcf.gz | bgzip -c > ${advancedFilteredDir}/${chr}_step2_depth_quality_filtered.vcf.gz
tabix -p vcf ${advancedFilteredDir}/${chr}_step2_depth_quality_filtered.vcf.gz

step2_count=$(bcftools view -H ${advancedFilteredDir}/${chr}_step2_depth_quality_filtered.vcf.gz | wc -l)
echo "SNP count after depth and quality filtering for ${chr}: $step2_count"

# Step 3: Filter SNPs near INDELs, multiallelic SNPs, INDELs, and monomorphic sites
bcftools filter -e 'AC==0 || AC==AN' --SnpGap 10 ${advancedFilteredDir}/${chr}_step2_depth_quality_filtered.vcf.gz | \
bcftools view -m2 -M2 -v snps -O v | bgzip -c > ${advancedFilteredDir}/${chr}_final_filtered.vcf.gz
tabix -p vcf ${advancedFilteredDir}/${chr}_final_filtered.vcf.gz

final_count=$(bcftools view -H ${advancedFilteredDir}/${chr}_final_filtered.vcf.gz | wc -l)
echo "Final SNP count for ${chr}: $final_count"

echo "Advanced filtering complete for chromosome ${chr}. Final filtered VCF: ${advancedFilteredDir}/${chr}_final_filtered.vcf.gz"

# Output SNP counts to a file
echo "${chr},${initial_count},${step1_count},${step2_count},${final_count}" >> ${advancedFilteredDir}/snp_counts.csv