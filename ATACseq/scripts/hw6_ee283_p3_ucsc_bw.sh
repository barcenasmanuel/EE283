#!/bin/bash

# Load required modules
module load samtools/1.15.1
module load bedtools2/2.30.0
module load ucsc-tools/v429

# Set directories and files
INPUT_DIR="/pub/mbarcen1/EE283/ATACseq/output/ucsc/sub_processed"
OUTPUT_DIR="/pub/mbarcen1/EE283/ATACseq/output/ucsc/custom_bigwig"
REF="/pub/mbarcen1/EE283/ref/dm6.fa.masked"
CHROM_SIZES="/pub/mbarcen1/EE283/ref/dm6.chrom.sizes"

# Set the condition (tissue type)
CONDITION="ED"  # Change this to "WD" for Wing Disc

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Combine BAM files for the condition
echo "Combining BAM files for $CONDITION..."
samtools merge -o $OUTPUT_DIR/${CONDITION}_combined.bam \
    $INPUT_DIR/A4_${CONDITION}_*.chrX.bam \
    $INPUT_DIR/A5_${CONDITION}_*.chrX.bam \
    $INPUT_DIR/A6_${CONDITION}_*.chrX.bam \
    $INPUT_DIR/A7_${CONDITION}_*.chrX.bam

# Sort and index the combined BAM file
samtools sort -o $OUTPUT_DIR/${CONDITION}_combined.sorted.bam $OUTPUT_DIR/${CONDITION}_combined.bam
samtools index $OUTPUT_DIR/${CONDITION}_combined.sorted.bam

# Count the total number of reads and create a scaling constant
echo "Calculating scaling factor..."
Nreads=$(samtools view -c -q 30 -F 4 $OUTPUT_DIR/${CONDITION}_combined.sorted.bam)
Scale=$(echo "1.0/($Nreads/1000000)" | bc -l)

# Calculate coverage and create bedGraph
echo "Calculating coverage..."
bedtools genomecov -ibam $OUTPUT_DIR/${CONDITION}_combined.sorted.bam -g $REF -bg -scale $Scale > $OUTPUT_DIR/${CONDITION}_coverage.bedgraph

# Convert bedGraph to bigWig
echo "Converting to bigWig..."
bedGraphToBigWig $OUTPUT_DIR/${CONDITION}_coverage.bedgraph $CHROM_SIZES $OUTPUT_DIR/${CONDITION}_custom.bw

echo "Custom bigWig file created: $OUTPUT_DIR/${CONDITION}_custom.bw"

# Clean up intermediate files
rm $OUTPUT_DIR/${CONDITION}_combined.bam $OUTPUT_DIR/${CONDITION}_combined.sorted.bam $OUTPUT_DIR/${CONDITION}_combined.sorted.bam.bai $OUTPUT_DIR/${CONDITION}_coverage.bedgraph

echo "Process completed."