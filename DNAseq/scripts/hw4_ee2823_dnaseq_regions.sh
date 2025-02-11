#!/bin/bash
#SBATCH --job-name=dna_extract_region_2    ## Name of the job
#SBATCH -A CLASS-ECOEVO283           ## Account to charge
#SBATCH -p standard                  ## Partition/queue name
#SBATCH --cpus-per-task=1            ## Number of cores the job needs
#SBATCH --mem=4G                     ## Memory required per node
#SBATCH --time=01:00:00              ## Time limit
#SBATCH --output=logs/dna_extract_region_2.%A.out
#SBATCH --error=logs/dna_extract_region_2.%A.err

module load samtools/1.10

# Define paths
outputDir="/pub/$USER/EE283/DNAseq/output"
region="X:1880000-2000000"  # Updated to match the chromosome naming in the BAM file
mapq_threshold=30

# Define the prefixes for the strains "A4" and "A5"
prefixes=("ADL06_1" "ADL06_2" "ADL06_3" "ADL09_1" "ADL09_2" "ADL09_3")

# Ensure output directory exists
mkdir -p $outputDir/regions

# Remove any previous output files
for prefix in "${prefixes[@]}"; do
  rm -f $outputDir/regions/${prefix}_X_1880000_2000000.bam
  rm -f $outputDir/regions/${prefix}_X_1880000_2000000.bam.bai
done

for prefix in "${prefixes[@]}"; do
  # Define input and output file names
  input_bam="${outputDir}/${prefix}.sorted.bam"
  output_bam="${outputDir}/regions/${prefix}_X_1880000_2000000.bam"

  # Check if the input BAM file exists
  if [[ ! -f $input_bam ]]; then
    echo "Warning: Input BAM file $input_bam not found. Skipping."
    continue
  fi

  # Extract the region with the specified mapq threshold
  samtools view -b -q $mapq_threshold $input_bam $region > $output_bam

  # Check if the output BAM file is not empty
  if [[ $(samtools view $output_bam | wc -l) -eq 0 ]]; then
    echo "Warning: No reads extracted for $prefix in the specified region with mapq > $mapq_threshold."
    continue
  fi

  # Index the output BAM file
  samtools index $output_bam

  echo "Processed ${input_bam} into ${output_bam}"
done