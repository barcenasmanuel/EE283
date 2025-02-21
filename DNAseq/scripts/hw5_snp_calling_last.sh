#!/bin/bash
#SBATCH --job-name=snp_calling
#SBATCH --account=CLASS-ECOEVO283
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --array=1-7
#SBATCH --output=logs/snp_calling_last_%A_%a.out
#SBATCH --error=logs/snp_calling_last_%A_%a.err

# Load required modules
module load java/1.8.0
module load gatk/4.2.6.1
module load samtools/1.15.1

# Define paths
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
inputDir="/pub/$USER/EE283/DNAseq/output"
outputDir="/pub/$USER/EE283/DNAseq/output/snp_results"
mkdir -p $outputDir

# Create chromosome names file if it doesn't exist
chromNamesFile="${outputDir}/chrome.names.txt"
if [[ ! -f $chromNamesFile ]]; then
    cat $ref | grep ">" | cut -f1 -d" " | tr -d ">" | head -n 7 > $chromNamesFile
fi

# Read sample prefixes from file
sampleFile="/pub/$USER/EE283/DNAseq/scripts/prefixes_DNAseq.txt"
samples=$(cat $sampleFile)

# Process each sample
for sample in $samples; do
    # Call SNPs with HaplotypeCaller
    if [[ ! -f ${outputDir}/${sample}.g.vcf.gz ]]; then
        /opt/apps/gatk/4.2.6.1/gatk HaplotypeCaller \
            -R $ref \
            -I ${inputDir}/${sample}.dedup.bam \
            --minimum-mapping-quality 30 \
            -ERC GVCF \
            -O ${outputDir}/${sample}.g.vcf.gz
    fi
done

# Combine GVCFs
if [[ ! -f ${outputDir}/allsample.g.vcf.gz ]]; then
    /opt/apps/gatk/4.2.6.1/gatk CombineGVCFs \
        -R $ref \
        $(printf -- '-V %s ' ${outputDir}/*.g.vcf.gz) \
        -O ${outputDir}/allsample.g.vcf.gz
fi

# Call SNPs by chromosome using array job
mychr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $chromNamesFile)

/opt/apps/gatk/4.2.6.1/gatk --java-options "-Xmx6g" GenotypeGVCFs \
    -R $ref \
    -V ${outputDir}/allsample.g.vcf.gz \
    --intervals $mychr \
    -stand-call-conf 5 \
    -O ${outputDir}/result.${mychr}.vcf.gz