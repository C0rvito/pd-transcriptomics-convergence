library(DESeq2)
library(ggplot2)
library(ggrepel)

# 1. Preparar os dados (assumindo que res_df agora tem a coluna gene_symbol)
# Se você carregou de um arquivo: 

setwd("/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/de")
res_df <- read.csv("4_SNCA_Triplication_counts_with_symbols.csv")
print(colnames(res_df))
# 2. Criar as categorias de expressão (conforme fizemos antes)
res_df$diff_expression <- "Não Significativo"
res_df$diff_expression[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated"
res_df$diff_expression[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down-regulated"

# 3. Criar uma coluna de label apenas para os genes que queremos mostrar
# Aqui selecionamos os top 10 genes com menor p-valor que são significativos
res_df$label <- ""
top_genes <- res_df[res_df$diff_expression != "Não Significativo", ]
top_genes <- top_genes[order(top_genes$padj), ]
top_labels <- head(top_genes$gene_id, 20) # Ajuste o número de labels aqui

res_df$label[rownames(res_df) %in% top_labels] <- as.character(res_df$gene_symbol[rownames(res_df) %in% top_labels]

# 4. Gerar o Volcano Plot com Labels
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
  geom_text_repel(aes(label = label), 
                  size = 3, 
                  max.overlaps = 20, 
                  box.padding = 0.5,
                  show.legend = FALSE) +
  labs(title = "Volcano Plot: SNCA Triplication vs Control",
       subtitle = "Top genes diferencialmente expressos com labels",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Status")

# 5. Salvar o gráfico
ggsave("Volcano_Plot_Annotated.png", p, width = 10, height = 8)

