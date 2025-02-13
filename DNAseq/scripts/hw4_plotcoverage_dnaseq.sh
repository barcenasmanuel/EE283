#!/bin/bash

# Activate the deeptools environment
mamba activate deeptools_ee283

# Define directories
inputDir="/pub/$USER/EE283/DNAseq/output/regions"
outputDir="/pub/$USER/EE283/DNAseq/output/coverage_plots"

# Ensure output directory exists
mkdir -p $outputDir

# Define prefixes
prefixes=("ADL06_1" "ADL06_2" "ADL06_3" "ADL09_1" "ADL09_2" "ADL09_3")

# Generate normalized coverage files using bamCoverage in bedGraph format
for prefix in "${prefixes[@]}"; do
    bamCoverage -b "${inputDir}/${prefix}_X_1880000_2000000.bam" \
                -o "${outputDir}/${prefix}_rpkm.bedgraph" \
                --outFileFormat bedgraph \
                --normalizeUsing RPKM
done

# Convert bedGraph to BED format if needed
# This step is optional and depends on the specific requirements of your downstream analysis
for prefix in "${prefixes[@]}"; do
    awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4}' "${outputDir}/${prefix}_rpkm.bedgraph" > "${outputDir}/${prefix}_rpkm.bed"
done

# Compute matrix for plotting using bedGraph files
computeMatrix reference-point -S "${outputDir}/ADL06_1_rpkm.bedgraph" "${outputDir}/ADL06_2_rpkm.bedgraph" "${outputDir}/ADL06_3_rpkm.bedgraph" \
                               "${outputDir}/ADL09_1_rpkm.bedgraph" "${outputDir}/ADL09_2_rpkm.bedgraph" "${outputDir}/ADL09_3_rpkm.bedgraph" \
                               -R "X:1880000:2000000" \
                               -a 5000 -b 5000 \
                               -o "${outputDir}/coverage_matrix.gz"

# Plot coverage
plotProfile -m "${outputDir}/coverage_matrix.gz" \
            --outFileName "${outputDir}/coverage_profile.png" \
            --regionsLabel "X:1,880,000-2,000,000" \
            --samplesLabel "A4_1" "A4_2" "A4_3" "A5_1" "A5_2" "A5_3" \
            --plotTitle "Coverage Profile for A4 and A5"