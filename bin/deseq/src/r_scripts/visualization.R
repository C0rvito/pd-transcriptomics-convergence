# >>>>>>>>>>>>>>>>>>>>>>>
# Arquivo: src/visualization.R
# Descrição: Gráficos de resultados de DE 
#   - Gráfico Volcano (LFC encolhido pelo apeglm)
#   - Gráfico MA
# <<<<<<<<<<<<<<<<<<<<<<<

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

source("bin/deseq/src/utils/utils.R")

# >>>>>>>>>>>>>>>>>>>>>>>
#' Gera um Gráfico Volcano a partir dos resultados encolhidos do DESeq2
#'
#' @param res_shrunken  data.frame com colunas: gene_symbol, log2FoldChange, padj.
#' @param target_name   Rótulo de texto para títulos e nomes de arquivos.
#' @param lfc_threshold Linhas tracejadas verticais (padrão 1).
#' @param padj_threshold Linha tracejada horizontal (padrão 0.05).
#' @param n_labels      Máximo de genes principais para rotular por direção (padrão 15).
#' @param plots_dir     Diretório de saída (padrão "plots").
#'
#' @return objeto ggplot (invisivelmente)
# <<<<<<<<<<<<<<<<<<<<<<<
generate_volcano_plot <- function(
    res_shrunken,
    target_name,
    lfc_threshold  = 1,
    padj_threshold = 0.05,
    n_labels       = 15,
    plots_dir      = "plots"
) {
  log_info("Gerando Gráfico Volcano: ", target_name)
  setup_dirs(plots_dir)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # Classificar genes
  # <<<<<<<<<<<<<<<<<<<<<<<
  res_df <- res_shrunken %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      Expressao = case_when(
        padj < padj_threshold & log2FoldChange >=  lfc_threshold ~ "Aumentado",
        padj < padj_threshold & log2FoldChange <= -lfc_threshold ~ "Reduzido",
        TRUE ~ "Não significativo"
      ),
      Expressao = factor(Expressao,
                          levels = c("Aumentado", "Reduzido", "Não significativo"))
    )

  # >>>>>>>>>>>>>>>>>>>>>>>
  # Selecionar genes para rotular: top n por padj dentro de cada grupo significativo
  # <<<<<<<<<<<<<<<<<<<<<<<
  genes_to_label <- res_df %>%
    dplyr::filter(Expressao != "Não significativo") %>%
    dplyr::group_by(Expressao) %>%
    dplyr::slice_min(order_by = padj, n = n_labels, with_ties = FALSE) %>%
    dplyr::ungroup()

  n_up   <- sum(res_df$Expressao == "Aumentado")
  n_down <- sum(res_df$Expressao == "Reduzido")
  subtitle <- paste0(
    "Mutant vs Control iPSCs  |  Aumentado: ", n_up,
    "  Reduzido: ", n_down,
    "  (|LFC| ≥ ", lfc_threshold, ", FDR < ", padj_threshold, ")"
  )

  p <- ggplot(res_df,
              aes(x = log2FoldChange, y = -log10(padj), color = Expressao)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_color_manual(
      values = c(
        "Aumentado"   = "firebrick",
        "Reduzido" = "steelblue",
        "Não significativo" = "grey60"
      )
    ) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", color = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = -log10(padj_threshold),
               linetype = "dashed", color = "grey30", linewidth = 0.5) +
    geom_text_repel(
      data          = genes_to_label,
      aes(label     = gene_symbol),
      size          = 3.2,
      max.overlaps  = 25,
      segment.color = "grey50",
      segment.size  = 0.3,
      show.legend   = FALSE
    ) +
    labs(
      title    = paste("Gráfico Volcano —", target_name),
      subtitle = subtitle,
      x        = "Log₂ Fold Change (encolhimento apeglm)",
      y        = "-Log₁₀ (P-valor Ajustado)",
      color    = "Expressão"
    ) +
    theme_pipeline() +
    guides(color = guide_legend(override.aes = list(size = 3)))

  out <- file.path(plots_dir, paste0("VolcanoPlot_", target_name, ".png"))
  ggsave(out, plot = p, width = 9, height = 6, dpi = 300)
  log_info("Gráfico Volcano salvo → ", out)

  invisible(p)
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Gera um Gráfico MA a partir dos resultados encolhidos do DESeq2
#'
#' @param res_shrunken  data.frame com colunas: gene_symbol, log2FoldChange,
#'                      padj, baseMean.
#' @param target_name   Rótulo de texto.
#' @param padj_threshold Limite de FDR para destaque (padrão 0.05).
#' @param n_labels       Máximo de genes para rotular (padrão 10).
#' @param plots_dir      Diretório de saída (padrão "plots").
#'
#' @return objeto ggplot (invisivelmente)
# <<<<<<<<<<<<<<<<<<<<<<<
generate_ma_plot <- function(
    res_shrunken,
    target_name,
    padj_threshold = 0.05,
    n_labels       = 10,
    plots_dir      = "plots"
) {
  log_info("Gerando Gráfico MA: ", target_name)
  setup_dirs(plots_dir)

  res_df <- res_shrunken %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange), baseMean > 0) %>%
    dplyr::mutate(Significativo = padj < padj_threshold & !is.na(padj))

  # >>>>>>>>>>>>>>>>>>>>>>>
  # Selecionar genes principais para rotular
  # <<<<<<<<<<<<<<<<<<<<<<<
  genes_to_label <- res_df %>%
    dplyr::filter(Significativo) %>%
    dplyr::slice_min(order_by = padj, n = n_labels, with_ties = FALSE)

  p <- ggplot(res_df, aes(x = log10(baseMean + 1), y = log2FoldChange,
                           color = Significativo)) +
    geom_point(alpha = 0.5, size = 1.2) +
    scale_color_manual(
      values = c("TRUE" = "firebrick", "FALSE" = "grey60"),
      labels = c("TRUE" = paste0("FDR < ", padj_threshold), "FALSE" = "Não significativo")
    ) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey30", linewidth = 0.5) +
    geom_text_repel(
      data         = genes_to_label,
      aes(label    = gene_symbol),
      size         = 3,
      max.overlaps = 20,
      show.legend  = FALSE
    ) +
    labs(
      title    = paste("Gráfico MA —", target_name),
      subtitle = "Mutant vs Control iPSCs (encolhimento de LFC apeglm)",
      x        = "Log₁₀ (Média das Contagens Normalizadas + 1)",
      y        = "Log₂ Fold Change",
      color    = NULL
    ) +
    theme_pipeline()

  out <- file.path(plots_dir, paste0("MAPlot_", target_name, ".png"))
  ggsave(out, plot = p, width = 8, height = 5, dpi = 300)
  log_info("Gráfico MA salvo → ", out)

  invisible(p)
}
