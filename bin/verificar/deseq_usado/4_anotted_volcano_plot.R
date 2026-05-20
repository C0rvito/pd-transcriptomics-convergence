library(DESeq2)
library(ggplot2)
library(ggrepel)

# 1. Preparar os dados

setwd("/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/de")
res_df <- read.csv("4_SNCA_Triplication_counts_with_symbols.csv")

# 2. Manter apenas linhas utilizaveis para o volcano plot
res_df <- subset(
  res_df,
  !is.na(log2FoldChange) & !is.na(pvalue) & !is.na(padj) & gene_symbol != ""
)

# 3. Criar as categorias de expressão
res_df$diff_expression <- "Não Significativo"
res_df$diff_expression[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated"
res_df$diff_expression[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down-regulated"

# 4. Selecionar apenas os genes mais diferencialmente expressos para rotular
top_n_labels <- 15
top_genes <- subset(res_df, diff_expression != "Não Significativo")
top_genes <- top_genes[order(-abs(top_genes$log2FoldChange), top_genes$padj), ]
top_genes <- head(top_genes, top_n_labels)

res_df$label <- ""
res_df$label[res_df$gene_id %in% top_genes$gene_id] <- as.character(
  res_df$gene_symbol[res_df$gene_id %in% top_genes$gene_id]
)

# 5. Gerar o Volcano Plot com Labels
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = diff_expression)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_minimal() +
  scale_color_manual(values = c("Down-regulated" = "blue", 
                                "Up-regulated" = "red", 
                                "Não Significativo" = "grey")) +
  # Linhas de referência
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  # Adicionar os nomes dos genes
  geom_text_repel(
                  data = subset(res_df, label != ""),
                  aes(label = label), 
                  size = 3, 
                  max.overlaps = 20, 
                  box.padding = 0.5,
                  show.legend = FALSE) +
  labs(title = "Volcano Plot: SNCA Triplication vs Control",
       subtitle = paste("Top", top_n_labels, "genes mais diferencialmente expressos com labels"),
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Status")

# 6. Salvar o gráfico
ggsave("Volcano_Plot_Annotated.png", p, width = 10, height = 8)

