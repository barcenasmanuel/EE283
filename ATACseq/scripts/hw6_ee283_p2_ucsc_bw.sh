#!/bin/bash

# Load required modules
module load ucsc-tools/v429
module load bedtools2/2.29.2

# Set directories
REF_DIR="/pub/mbarcen1/EE283/ref"
INPUT_DIR="/pub/mbarcen1/EE283/ATACseq/output/ucsc/peak_calling"
OUTPUT_DIR="/pub/mbarcen1/EE283/ATACseq/output/ucsc/bigwig"

# Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Check if chromosome sizes file exists, if not create it
if [ ! -f "$REF_DIR/dm6.chrom.sizes" ]; then
    echo "Creating chromosome sizes file..."
    cut -f1,2 $REF_DIR/dm6.fa.masked.fai > $REF_DIR/dm6.chrom.sizes
fi

# Convert BED to bedGraph
echo "Converting BED to bedGraph..."
bedtools genomecov -i $INPUT_DIR/ED_combined.tn5.bed -g $REF_DIR/dm6.chrom.sizes -bg > $OUTPUT_DIR/ED_treat_pileup.bdg

# Process the bedGraph file
echo "Processing bedGraph file..."
LC_COLLATE=C sort -k1,1 -k2,2n $OUTPUT_DIR/ED_treat_pileup.bdg > $OUTPUT_DIR/ED_treat_pileup.sorted.bdg

# Trim extended fragments
echo "Trimming extended fragments..."
awk 'NR==FNR {chr_size[$1]=$2; next} $3 <= chr_size[$1]' $REF_DIR/dm6.chrom.sizes $OUTPUT_DIR/ED_treat_pileup.sorted.bdg > $OUTPUT_DIR/ED_treat_pileup.safeends.bdg

# Convert to bigWig
echo "Converting to bigWig..."
bedGraphToBigWig $OUTPUT_DIR/ED_treat_pileup.safeends.bdg $REF_DIR/dm6.chrom.sizes $OUTPUT_DIR/ED_broad_peaks.bw

echo "BigWig file created: $OUTPUT_DIR/ED_broad_peaks.bw"