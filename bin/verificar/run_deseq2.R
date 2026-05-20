library(DESeq2)
library(ggplot2)

# Set working directory
setwd("/Users/mateuslisboa/Desktop/Data/LRRK2")

# 1. Load data
counts <- read.csv("LRRK2_mRNA_counts_matrix.csv", row.names = 1)
metadata <- read.csv("metadata.csv", row.names = 1)

# Ensure sample names match
counts <- counts[, rownames(metadata)]

run_dea <- function(subset_metadata, comparison_name) {
  print(paste("Running DEA for:", comparison_name))
  
  # Subset counts
  subset_counts <- counts[, rownames(subset_metadata)]
  
  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(countData = subset_counts,
                                colData = subset_metadata,
                                design = ~ condition)
  
  # Set reference level
  dds$condition <- relevel(dds$condition, ref = "Control")
  
  # Run DESeq
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds)
  res_df <- as.data.frame(res)
  
  # Add gene symbols
  # (Since we have the matrix with symbols, let's join it or just use the IDs)
  # For now, we'll just save the IDs.
  
  # Save results
  write.csv(res_df, paste0("DESeq2_Results_", comparison_name, ".csv"))
  
  # Simple Volcano Plot
  res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Significant", "Not Significant")
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point(alpha = 0.4) +
    theme_minimal() +
    scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
    labs(title = paste("Volcano Plot:", comparison_name),
         x = "Log2 Fold Change",
         y = "-Log10 P-value")
  
  ggsave(paste0("Volcano_Plot_", comparison_name, ".png"), p, width = 8, height = 6)
  
  return(dds)
}

# 2. Comparison 1: Non-Isogenic
meta_non_iso <- metadata[metadata$comparison == "Non-Isogenic", ]
dds_non_iso <- run_dea(meta_non_iso, "Non_Isogenic")

# 3. Comparison 2: Isogenic
meta_iso <- metadata[metadata$comparison == "Isogenic", ]
dds_iso <- run_dea(meta_iso, "Isogenic")

print("DESeq2 Analysis Complete.")
