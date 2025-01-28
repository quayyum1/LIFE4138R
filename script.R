#Summary Statistics

# Load necessary libraries
library(dplyr)

# Load and read TSV files
AvsF <- read.delim("A_vs_F.deseq2.results.tsv")
AvsC <- read.delim("A_vs_C.deseq2.results.tsv")

# Give thresholds for significance
padj_threshold <- 0.05
logfold_threshold <- 1.0

# Function for comparison anaylsis
analyse_comparison <- function(df, comparison_name) {
  upregulated <- nrow(df[df$padj < padj_threshold & df$log2FoldChange > logfold_threshold, ]) # Set the requirements to class upregulation
  downregulated <- nrow(df[df$padj < padj_threshold & df$log2FoldChange < -logfold_threshold, ]) # Set the requirements to class downregulation
  
  padj_summary <- summary(df$padj) # Summary function for adjusted p Value
  logfc_summary <- summary(df$log2FoldChange)
  
  cat("\nComparison:", comparison_name, "\n")
  cat("Significantly Upregulated Genes:", upregulated, "\n")
  cat("Significantly Downregulated Genes:", downregulated, "\n")
  cat("Adjusted P-Value (padj) Summary:\n")
  print(padj_summary)
  cat("Log2 Fold Change Summary:\n")
  print(logfc_summary)
}

# Analyse both comparisons
summary_AvsF <- analyse_comparison(AvsF, "A_vs_F")
summary_AvsC <- analyse_comparison(AvsC, "A_vs_C")

#===================================================================================================================================================
# Plots
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)

volcano_plot <- function(df, comparison_name, padj_threshold = 0.05, logfold_threshold = 1) {
  # Filter out rows with NA values or invalid values
  df <- df %>% filter(!is.na(padj) & !is.na(log2FoldChange) & is.finite(padj) & is.finite(log2FoldChange))
  
  # Create the volcano plot
  ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = padj < padj_threshold & abs(log2FoldChange) > logfold_threshold)) +
    geom_point(alpha = 0.6) + # Add points with transparency
    scale_color_manual(values = c("black", "red")) + # Highlight significant points in red
    labs(title = paste("Volcano Plot -", comparison_name), x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value") +
    theme_minimal() # Use a clean minimal theme
}

# MA plot function
# Displays the relationship between mean expression (baseMean) and log2 fold change
ma_plot <- function(df, comparison_name) {
  ggplot(df, aes(x = baseMean, y = log2FoldChange, color = padj < padj_threshold)) +
    geom_point(alpha = 0.6) + # Add points with transparency
    scale_x_log10() + # Log-transform the x-axis to better display large ranges
    labs(title = paste("MA Plot -", comparison_name), x = "Mean Expression (log10)", y = "Log2 Fold Change") +
    theme_minimal() # Use a clean minimal theme
}

# Histogram of p-values
pvalue_histogram <- function(df, comparison_name) {
  ggplot(df, aes(x = padj)) +
    geom_histogram(binwidth = 0.01, fill = "blue", alpha = 0.7) + # Create bins of width 0.01
    labs(title = paste("P-Value Histogram -", comparison_name), x = "Adjusted P-Value", y = "Frequency") +
    theme_minimal()
}

# Heatmap of top expressed genes with smallest adjusted
heatmap_plot <- function(df, comparison_name) {
  top_genes <- df %>% 
    filter(padj < padj_threshold) %>% # Filter genes based on significance threshold
    arrange(padj) %>% # Sort by adjusted p-value
    head(50) %>% # Select the top 50 genes
    select(log2FoldChange) # Extract the log2 fold change column
  mat <- as.matrix(top_genes) # Convert the data to a matrix
  rownames(mat) <- paste("Gene", 1:nrow(mat)) # Assign row names to identify genes
  pheatmap(mat, cluster_rows = TRUE, cluster_cols = FALSE, # Cluster rows but not columns
           main = paste("Heatmap of Top Genes -", comparison_name), 
           color = colorRampPalette(c("blue", "white", "red"))(50)) # Use a blue-white-red color gradient
}

# Generate plots for A_vs_F
# Call each plotting function for the first dataset
volcano_plot(AvsF, "A_vs_F")
print(ma_plot(AvsF, "A_vs_F"))
pvalue_histogram(AvsF, "A_vs_F")
heatmap_plot(AvsF, "A_vs_F")

# Generate plots for A_vs_C
# Call each plotting function for the second dataset
volcano_plot(AvsC, "A_vs_C")
ma_plot(AvsC, "A_vs_C")
pvalue_histogram(AvsC, "A_vs_C")
heatmap_plot(AvsC, "A_vs_C")

ma_AvsF <- ma_plot(AvsF, "A_vs_F")
hist_AvsF <- pvalue_histogram(AvsF, "A_vs_F")
heatmap_AvsF <- heatmap_plot(AvsF, "A_vs_F")

ma_AvsC <- ma_plot(AvsC, "A_vs_C")
hist_AvsC <- pvalue_histogram(AvsC, "A_vs_C")
heatmap_AvsC <- heatmap_plot(AvsC, "A_vs_C")

