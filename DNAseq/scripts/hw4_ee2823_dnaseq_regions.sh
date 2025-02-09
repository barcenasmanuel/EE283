#!/bin/bash
#SBATCH --job-name=extract_region    ## Name of the job
#SBATCH -A CLASS-ECOEVO283           ## Account to charge
#SBATCH -p standard                  ## Partition/queue name
#SBATCH --cpus-per-task=1            ## Number of cores the job needs
#SBATCH --mem=4G                     ## Memory required per node
#SBATCH --time=01:00:00              ## Time limit
#SBATCH --output=logs/extract_region.%A.out
#SBATCH --error=logs/extract_region.%A.err

module load samtools/1.10

# Define paths
outputDir="/pub/$USER/EE283/DNAseq/output"
region="chrX:1880000-2000000"
mapq_threshold=30

# Define the prefixes for the strains "A4" and "A5"
prefixes=("ADL06_1" "ADL06_2" "ADL06_3" "ADL09_1" "ADL09_2" "ADL09_3")

# Ensure output directory exists
mkdir -p $outputDir/regions

# Remove any previous output files
for prefix in "${prefixes[@]}"; do
  rm -f $outputDir/regions/${prefix}_chrX_1880000_2000000.bam
  rm -f $outputDir/regions/${prefix}_chrX_1880000_2000000.bam.bai
done

for prefix in "${prefixes[@]}"; do
  # Define input and output file names
  input_bam="${outputDir}/${prefix}.sorted.bam"
  output_bam="${outputDir}/regions/${prefix}_chrX_1880000_2000000.bam"

  # Check if the input BAM file exists
  if [[ ! -f $input_bam ]]; then
    echo "Warning: Input BAM file $input_bam not found. Skipping."
    continue
  fi

  # Extract the region with the specified mapq threshold
  samtools view -b -q $mapq_threshold $input_bam $region > $output_bam

  # Index the output BAM file
  samtools index $output_bam

  echo "Processed ${input_bam} into ${output_bam}"
done