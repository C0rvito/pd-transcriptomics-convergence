library(DESeq2)
library(ggplot2)

# 1. Configurar diretório de trabalho
setwd("/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/de")

# 2. Carregar dados
counts <- read.csv("matriz_contagem_final.csv", row.names = 1)
metadata <- read.csv("~/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/samples/deseq_metadata.csv", row.names = 1)

# 3. Sincronizar amostras
# Garante que a matriz de contagem tenha as mesmas colunas que as linhas do metadata
counts <- counts[, rownames(metadata)]

# 4. Criar o objeto DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ source_name)

# 5. Definir o nível de referência
dds$source_name <- relevel(dds$source_name, ref = "Control_NSC")

# 6. Executar a análise de Expressão Diferencial
dds <- DESeq(dds)

# 7. Extrair resultados: Triplicação SNCA vs Controle (ambos em No-Seeded-Control)
res <- results(dds, contrast=c("source_name", "Triplication SNCA_NSC", "Control_NSC"))
res_df <- as.data.frame(res)

# 8. Salvar os resultados em CSV
write.csv(res_df, "DESeq2_SNCA_Triplication_vs_Control_Basal.csv")

# 9. Categorizar os genes para o Volcano Plot
res_df$diff_expression <- "Não Significativo"

# Genes Up-regulated (Log2FC > 1 e p-adj < 0.05)
res_df$diff_expression[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated"

# Genes Down-regulated (Log2FC < -1 e p-adj < 0.05)
res_df$diff_expression[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down-regulated"

# 10. Criar o Volcano Plot colorido
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = diff_expression)) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_minimal() +
  # Definindo cores: Azul para Down e Vermelho para Up
  scale_color_manual(values = c("Down-regulated" = "blue", 
                                "Up-regulated" = "red", 
                                "Não Significativo" = "grey")) +
  # Linhas de referência para os cortes (thresholds)
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  labs(title = "Volcano Plot: Triplicação SNCA vs Controle",
       subtitle = "Genes Up e Down regulated (p-adj < 0.05 | |Log2FC| > 1)",
       x = "Log2 Fold Change",
       y = "-Log10 P-value",
       color = "Status")

# 11. Salvar o gráfico final
ggsave("Volcano_Plot_SNCA.png", p, width = 8, height = 6)

print("Análise concluída. Resultados e gráfico salvos com sucesso!")

# 12 Filtrando Genes Diferencialmente Expressos
degs <- res_df[which(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1), ]
write.csv(degs, "lista_DEGs_significativos.csv")
