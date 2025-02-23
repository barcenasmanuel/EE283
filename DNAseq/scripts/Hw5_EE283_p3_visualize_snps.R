#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)

# Set the working directory
setwd("/pub/mbarcen1/EE283/DNAseq/output/snp_results/visualization_advanced")

# Read the 012 file
data <- read.table("X_1Mb_012.012", header=FALSE)

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
  labs(title = "Biallelic SNPs in first 1Mb of X chromosome (Advanced Filtered)",
       x = "SNPs", y = "Samples") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# Save the plot
ggsave("X_1Mb_snps_visualization_advanced.png", p, width = 12, height = 6, dpi = 300)

print("Visualization complete. Output image: X_1Mb_snps_visualization_advanced.png")