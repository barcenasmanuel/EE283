import matplotlib.pyplot as plt
import pandas as pd

# Define the path to the bedGraph files and the output directory
inputDir = '/pub/mbarcen1/EE283/DNAseq/output/coverage_plots'
outputDir = '/pub/mbarcen1/EE283/DNAseq/output/coverage_plots'

# Load the RPKM data for the samples you want to compare
adl06_1_coverage = pd.read_csv(f'{inputDir}/ADL06_1_rpkm.bedgraph', sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])
adl09_1_coverage = pd.read_csv(f'{inputDir}/ADL09_1_rpkm.bedgraph', sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])

# Focus on the region around 1,904,042
region_start = 1880000
region_end = 2000000

adl06_1_region = adl06_1_coverage[(adl06_1_coverage['start'] >= region_start) & (adl06_1_coverage['end'] <= region_end)]
adl09_1_region = adl09_1_coverage[(adl09_1_coverage['start'] >= region_start) & (adl09_1_coverage['end'] <= region_end)]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(adl06_1_region['start'], adl06_1_region['coverage'], label='ADL06_1', color='blue')
plt.plot(adl09_1_region['start'], adl09_1_region['coverage'], label='ADL09_1', color='red')
plt.axvline(x=1904042, color='green', linestyle='--', label='Position 1,904,042')
plt.xlabel('Genomic Position')
plt.ylabel('RPKM Coverage')
plt.title('Per Base Pair Read Coverage')
plt.legend()

# Save the plot to the output directory
output_file = f'{outputDir}/coverage_profile_ADL06_1_ADL09_1.png'
plt.savefig(output_file, dpi=300, bbox_inches='tight')
plt.close()

print(f"Plot saved to {output_file}")