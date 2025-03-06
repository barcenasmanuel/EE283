#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH -A CLASS-ECOEVO283
#SBATCH -p standard
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/featureCounts.out
#SBATCH --error=logs/featureCounts.err

# Load the required module
module load subread/2.0.3

# Set the paths
gtf="/pub/mbarcen1/EE283/ref/dmel-all-r6.13.gtf"
output_dir="/pub/mbarcen1/EE283/RNAseq/counts"

# Create output directory if it doesn't exist
mkdir -p $output_dir

# Create a space-delimited list of BAM files
myfile=$(cat /pub/mbarcen1/EE283/hw7/shortRNAseq.names.txt | tr "\n" " ")

# Run featureCounts
featureCounts -p -T 8 -t exon -g gene_id -Q 30 -F GTF -a $gtf -o $output_dir/fly_counts.txt $myfile

echo "featureCounts completed. Output saved to $output_dir/fly_counts.txt"