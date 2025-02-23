#!/usr/bin/env Rscript

# Load required libraries
library(ggplot2)

# Read the SNP count data
data <- read.table("/pub/mbarcen1/EE283/DNAseq/output/snp_results/advanced_filtered_snps/snp_counts.csv", 
                   header = FALSE, sep = ",",
                   col.names = c("Chromosome", "Initial", "Step1", "Step2", "Final"))

# Calculate total SNPs at each step
total_snps <- colSums(data[, -1])  # Exclude the Chromosome column
print("Total SNPs at each step:")
print(total_snps)

# Calculate SNPs lost at each step
snps_lost <- c(0, diff(total_snps) * -1)
names(snps_lost) <- names(total_snps)

# Create a data frame for plotting
plot_data <- data.frame(
  Step = factor(names(total_snps), levels = names(total_snps)),
  SNPs_Remaining = total_snps,
  SNPs_Lost = snps_lost
)

# Create the plot
p <- ggplot(plot_data, aes(x = Step)) +
  geom_col(aes(y = SNPs_Remaining), fill = "blue", alpha = 0.7) +
  geom_col(aes(y = SNPs_Lost), fill = "red", alpha = 0.7) +
  geom_text(aes(y = SNPs_Remaining, label = scales::comma(SNPs_Remaining)), 
            vjust = -0.5, color = "blue", size = 3) +
  geom_text(aes(y = SNPs_Lost, label = scales::comma(SNPs_Lost)), 
            vjust = 1.5, color = "red", size = 3) +
  theme_minimal() +
  labs(title = "SNP Counts and Losses After Each Filtering Step",
       x = "Filtering Step",
       y = "Number of SNPs") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(labels = scales::comma) +
  annotate("text", x = 1, y = max(total_snps), 
           label = paste("Initial SNPs:", format(total_snps[1], big.mark = ",")), 
           hjust = 0, vjust = 1, color = "black", fontface = "bold") +
  annotate("text", x = length(total_snps), y = total_snps[length(total_snps)] * 1.3, 
           label = paste("Final SNPs:", format(total_snps[length(total_snps)], big.mark = ",")), 
           hjust = 1, vjust = -0.5, color = "black", fontface = "bold")

# Save the plot
ggsave("/pub/mbarcen1/EE283/DNAseq/output/snp_results/advanced_filtered_snps/snp_filtering_loss_plot.png", 
       p, width = 12, height = 8, dpi = 300)

# Calculate and print the percentage of SNPs retained
percent_retained <- (total_snps["Final"] / total_snps["Initial"]) * 100

cat(sprintf("\nPercentage of SNPs retained: %.2f%%\n", percent_retained))

# Print the percentage change between each step
steps <- names(total_snps)
for (i in 2:length(steps)) {
  percent_change <- (total_snps[steps[i]] - total_snps[steps[i-1]]) / total_snps[steps[i-1]] * 100
  cat(sprintf("Percentage change from %s to %s: %.2f%%\n", steps[i-1], steps[i], percent_change))
}

cat("\nNote: The number of SNPs at each step is estimated from the VCF files.\n")
cat("The blue bars show the remaining SNPs, while the red bars show the SNPs lost at each step.\n")