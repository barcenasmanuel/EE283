#!/bin/bash

module load samtools/1.15.1
module load bedtools2/2.30.0
module load macs/2.2.7.1

# Set the input and output directories
INPUT_DIR="/pub/mbarcen1/EE283/ATACseq/output/ucsc/sub_processed"
OUTPUT_DIR="/pub/mbarcen1/EE283/ATACseq/output/ucsc/peak_calling"

# Set the condition name
CONDITION="ED"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

echo "Starting analysis for condition: $CONDITION"

# Step 1: Merge BAM files for the condition
echo "Merging BAM files..."
samtools merge -o $OUTPUT_DIR/${CONDITION}_combined.bam \
    $INPUT_DIR/A4_${CONDITION}_*.chrX.bam \
    $INPUT_DIR/A5_${CONDITION}_*.chrX.bam \
    $INPUT_DIR/A6_${CONDITION}_*.chrX.bam \
    $INPUT_DIR/A7_${CONDITION}_*.chrX.bam

# Step 2: Convert to BED with Tn5 offset correction
echo "Converting to BED and applying Tn5 offset correction..."
bedtools bamtobed -i $OUTPUT_DIR/${CONDITION}_combined.bam | \
awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' \
> $OUTPUT_DIR/${CONDITION}_combined.tn5.bed

# Step 3: Call peaks using MACS2
echo "Calling peaks with MACS2..."
macs2 callpeak -t $OUTPUT_DIR/${CONDITION}_combined.tn5.bed \
    -n ${CONDITION} \
    -f BED \
    -g dm \
    -q 0.01 \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    --call-summits \
    --keep-dup all \
    -B \
    --broad \
    --outdir $OUTPUT_DIR

echo "Peak calling completed for ${CONDITION}. Output files are in $OUTPUT_DIR"