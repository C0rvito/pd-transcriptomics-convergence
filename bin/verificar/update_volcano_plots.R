library(ggplot2)
library(ggrepel)
library(dplyr)

# Set working directory
setwd("/Users/mateuslisboa/Desktop/Data/LRRK2")

results_dir <- "outputs/differential_expression"
plots_dir <- "outputs/plots"

files_to_process <- c(
  "Isogenic" = "DESeq2_Results_Isogenic_Annotated.csv",
  "Non_Isogenic" = "DESeq2_Results_Non_Isogenic_Annotated.csv",
  "Combined" = "DESeq2_Results_Combined_Annotated.csv"
)

for (label in names(files_to_process)) {
  file_path <- file.path(results_dir, files_to_process[label])
  
  if (!file.exists(file_path)) {
    print(paste("File not found:", file_path))
    next
  }
  
  print(paste("Processing volcano plot for:", label))
  
  # Load data
  res_df <- read.csv(file_path)
  
  # Remove rows with NA in padj or pvalue
  res_df <- res_df %>% filter(!is.na(padj) & !is.na(pvalue))
  
  # Define significance and color
  res_df$diffexpressed <- "Not Significant"
  res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated"
  res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down-regulated"
  
  # Identify top 10 up and top 10 down
  top_up <- res_df %>% 
    filter(diffexpressed == "Up-regulated") %>% 
    arrange(desc(log2FoldChange)) %>% 
    head(10)
    
  top_down <- res_df %>% 
    filter(diffexpressed == "Down-regulated") %>% 
    arrange(log2FoldChange) %>% 
    head(10)
    
  # Combine labels
  label_df <- rbind(top_up, top_down)
  
  # Volcano Plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = diffexpressed)) +
    geom_point(alpha = 0.5, size = 1.5) +
    scale_color_manual(values = c("Down-regulated" = "blue", 
                                 "Not Significant" = "grey", 
                                 "Up-regulated" = "red")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
    theme_minimal() +
    labs(title = paste("Volcano Plot:", label),
         subtitle = "Labels: Top 10 Up and Top 10 Down Regulated Genes",
         x = "log2 Fold Change",
         y = "-log10 p-value",
         color = "Expression Status") +
    theme(legend.position = "bottom") +
    # Add labels
    geom_text_repel(data = label_df, aes(label = gene_symbol),
                    size = 3.5,
                    box.padding = 0.5,
                    point.padding = 0.3,
                    max.overlaps = Inf,
                    show.legend = FALSE)

  # Save plot
  output_plot <- file.path(plots_dir, paste0("Volcano_Plot_", label, "_Updated.png"))
  ggsave(output_plot, p, width = 10, height = 8, dpi = 300)
  print(paste("Saved:", output_plot))
}

print("All volcano plots updated.")
