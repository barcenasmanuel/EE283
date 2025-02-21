#!/bin/bash
#SBATCH --job-name=snp_calling_pipeline
#SBATCH --account=CLASS-ECOEVO283
#SBATCH --partition=standard
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/snp_calling.%A.out
#SBATCH --error=logs/snp_calling.%A.err

module load java/1.8.0
module load gatk/4.2.6.1
module load picard-tools/2.27.1
module load samtools/1.15.1

# Define paths
ref="/pub/$USER/EE283/ref/dmel-all-chromosome-r6.13.fasta"
inputDir="/pub/$USER/EE283/DNAseq/output"
outputDir="/pub/$USER/EE283/DNAseq/output/snp_results"
mkdir -p $outputDir

# Read genotype prefixes from file and extract unique genotypes
genotypeFile="/pub/$USER/EE283/DNAseq/scripts/prefixes_DNAseq.txt"
genotypes=$(cut -d'_' -f1 $genotypeFile | sort | uniq)

# Process each genotype
for genotype in $genotypes; do
    # Find and merge BAM files for the genotype
    bamFiles=$(ls ${inputDir}/${genotype}_*.sorted.bam)
    samtools merge -o ${outputDir}/${genotype}.bam $bamFiles
    samtools sort ${outputDir}/${genotype}.bam -o ${outputDir}/${genotype}.sort.bam

    # Add read groups
    java -jar /opt/apps/picard-tools/2.27.1/picard.jar AddOrReplaceReadGroups \
        I=${outputDir}/${genotype}.sort.bam \
        O=${outputDir}/${genotype}.RG.bam \
        SORT_ORDER=coordinate \
        RGPL=illumina RGPU=D109LACXX RGLB=Lib1 RGID=${genotype} RGSM=${genotype} \
        VALIDATION_STRINGENCY=LENIENT

    # Remove duplicates
    java -jar /opt/apps/picard-tools/2.27.1/picard.jar MarkDuplicates REMOVE_DUPLICATES=true \
        I=${outputDir}/${genotype}.RG.bam \
        O=${outputDir}/${genotype}.dedup.bam \
        M=${outputDir}/${genotype}_marked_dup_metrics.txt

    samtools index ${outputDir}/${genotype}.dedup.bam

    # Call SNPs with HaplotypeCaller
    /opt/apps/gatk/4.2.6.1/gatk HaplotypeCaller \
        -R $ref \
        -I ${outputDir}/${genotype}.dedup.bam \
        --minimum-mapping-quality 30 \
        -ERC GVCF \
        -O ${outputDir}/${genotype}.g.vcf.gz
done

# Combine GVCFs
/opt/apps/gatk/4.2.6.1/gatk CombineGVCFs \
    -R $ref \
    $(printf -- '-V %s ' ${outputDir}/*.g.vcf.gz) \
    -O ${outputDir}/allsample.g.vcf.gz

# Prepare chromosome names file
cat $ref | grep ">" | cut -f1 -d" " | tr -d ">" | head -n 7 > ${outputDir}/chrome.names.txt

# Call SNPs by chromosome using an array job
mychr=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${outputDir}/chrome.names.txt)

/opt/apps/gatk/4.2.6.1/gatk --java-options "-Xmx6g" GenotypeGVCFs \
    -R $ref \
    -V ${outputDir}/allsample.g.vcf.gz \
    --intervals $mychr \
    -stand-call-conf 5 \
    -O ${outputDir}/result.${mychr}.vcf.gz