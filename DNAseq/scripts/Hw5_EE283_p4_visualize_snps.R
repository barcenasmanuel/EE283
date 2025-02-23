#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)

# Set the working directory
setwd("/pub/mbarcen1/EE283/DNAseq/output/snp_results/visualization_task4_advanced")

# Read the filtered 012 file
data <- read.table("X_1Mb_012_filtered_task4.012", header=FALSE)

# Extract the genotype data (excluding the first two columns)
genotypes <- as.matrix(data[, 3:ncol(data)])

# Transpose the matrix
genotypes_t <- t(genotypes)

# Create a data frame for ggplot
df <- expand.grid(SNP = 1:nrow(genotypes_t), Sample = 1:ncol(genotypes_t))
df$Genotype <- as.vector(genotypes_t)

# Create the plot
p <- ggplot(df, aes(x = SNP, y = Sample, fill = factor(Genotype))) +
  geom_tile() +
  scale_fill_manual(values = c("grey", "green", "yellow", "red"),
                    labels = c("Missing", "Ref/Ref", "Ref/Alt", "Alt/Alt"),
                    name = "Genotype") +
  theme_minimal() +
  labs(title = "Filtered Biallelic SNPs in first 1Mb of X chromosome",
       subtitle = "SNPs with two strains 0/0 and two strains 1/1",
       x = "SNPs", y = "Samples") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save the plot
ggsave("X_1Mb_snps_visualization_task4_advanced.png", p, width = 12, height = 6, dpi = 300)

print("Task 4 Advanced: Visualization complete. Output image: X_1Mb_snps_visualization_task4_advanced.png")

# Print summary statistics
total_snps <- nrow(genotypes_t)
cat(sprintf("Total SNPs after filtering: %d\n", total_snps))

# Count occurrences of each genotype
genotype_counts <- table(df$Genotype)
cat("Genotype counts:\n")
print(genotype_counts)

# Calculate percentage of SNPs with 50% frequency (assuming 0 and 2 represent the two homozygous states)
snps_50_percent <- sum(rowSums(genotypes_t == 0) == 2 & rowSums(genotypes_t == 2) == 2)
cat(sprintf("SNPs with exactly two 0/0 and two 2/2: %d\n", snps_50_percent))

# Try to read the total number of SNPs from the VCF file
tryCatch({
  total_vcf_snps <- as.numeric(system("bcftools view -H X_1Mb_task4.vcf.gz | wc -l", intern = TRUE))
  cat(sprintf("Percentage of SNPs with 50%% frequency: %.2f%%\n", (snps_50_percent / total_vcf_snps) * 100))
}, error = function(e) {
  cat("Error reading total SNP count from VCF file. Make sure bcftools is available and the file exists.\n")
})