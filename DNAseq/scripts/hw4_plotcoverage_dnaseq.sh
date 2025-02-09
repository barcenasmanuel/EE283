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

# Generate normalized coverage files using bamCoverage
for prefix in "${prefixes[@]}"; do
    bamCoverage -b ${inputDir}/${prefix}_chrX_1880000_2000000.bam \
                -o ${outputDir}/${prefix}_rpkm.bw \
                --normalizeUsing RPKM
done

# Compute matrix for plotting
computeMatrix reference-point -S ${outputDir}/ADL06_1_rpkm.bw ${outputDir}/ADL06_2_rpkm.bw ${outputDir}/ADL06_3_rpkm.bw \
                               ${outputDir}/ADL09_1_rpkm.bw ${outputDir}/ADL09_2_rpkm.bw ${outputDir}/ADL09_3_rpkm.bw \
                               -R chrX:1880000:2000000 \
                               -a 5000 -b 5000 \
                               -o ${outputDir}/coverage_matrix.gz

# Plot coverage
plotProfile -m ${outputDir}/coverage_matrix.gz \
            --outFileName ${outputDir}/coverage_profile.png \
            --regionsLabel "chrX:1,880,000-2,000,000" \
            --samplesLabel "A4_1" "A4_2" "A4_3" "A5_1" "A5_2" "A5_3" \
            --plotTitle "Coverage Profile for A4 and A5"