# >>>>>>>>>>>>>>>>>>>>>>>
# Arquivo: src/gsea_pipeline.R
# Descrição: Gene Ontology (ORA) e Análise de Enriquecimento de Conjunto de Genes (GSEA)
#   - ORA via enrichGO em DEGs significativos (com proteções para listas pequenas)
#   - GSEA via gseGO ranqueado pela estatística Wald (mais poderoso que log2FC isolado)
#   - Simplificação de termos GO para reduzir redundância (similaridade semântica de Wang)
#   - Exporta resultados e gera dotplots, barplots, ridge plots
# <<<<<<<<<<<<<<<<<<<<<<<

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(ggridges)
  library(scales)
})

source("bin/deseq/src/utils/utils.R")

# >>>>>>>>>>>>>>>>>>>>>>>
# Mínimo de DEGs necessários para tentar ORA (abaixo disso, os resultados não são confiáveis)
# <<<<<<<<<<<<<<<<<<<<<<<
.MIN_GENES_ORA <- 5

# >>>>>>>>>>>>>>>>>>>>>>>
#' Executa a Análise de Sobre-representação de GO (ORA) e GSEA
#'
#' @param res_shrunken  data.frame (apeglm): gene_symbol, log2FoldChange, padj.
#'                      Usado para rótulos de volcano e lista de genes ORA.
#' @param res_raw       data.frame (resultados brutos do DESeq2): deve incluir `stat`
#'                      (estatística Wald). Usado para ranquear todos os genes para o GSEA.
#' @param degs          data.frame de DEGs significativos (pré-filtrados).
#' @param target_name   Rótulo de texto para títulos e nomes de arquivos.
#' @param padj_ora      Corte de FDR para GO ORA (padrão 0.05).
#' @param padj_gsea     Corte de FDR para GSEA (padrão 0.05).
#' @param min_gs_size   Tamanho mínimo do conjunto de genes para GSEA (padrão 15).
#' @param max_gs_size   Tamanho máximo do conjunto de genes para GSEA (padrão 500).
#' @param out_dir       Diretório de saída para resultados TSV (padrão "results").
#' @param plots_dir     Diretório de saída para gráficos (padrão "plots").
#'
#' @return Lista nomeada: ego (resultado ORA), gse (resultado GSEA), ou NULL se pulado.
# <<<<<<<<<<<<<<<<<<<<<<<
run_enrichment_analysis <- function(
    res_shrunken,
    res_raw,
    degs,
    target_name,
    padj_ora    = 0.05,
    padj_gsea   = 0.05,
    min_gs_size = 15,
    max_gs_size = 500,
    out_dir     = "results",
    plots_dir   = "plots"
) {
  log_info("=== Análise de Enriquecimento: ", target_name, " ===")
  setup_dirs(c(out_dir, plots_dir))

  results_list <- list(ego = NULL, gse = NULL)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 1. Análise de Sobre-representação de GO (ORA)
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("--- Passo 1: GO ORA ---")
  genes_sig <- degs %>%
    dplyr::filter(!is.na(gene_symbol)) %>%
    dplyr::pull(gene_symbol) %>%
    unique()

  log_info("DEGs disponíveis para ORA: ", length(genes_sig))

  if (length(genes_sig) < .MIN_GENES_ORA) {
    log_warn("Apenas ", length(genes_sig), " DEGs encontrados. ",
             "GO ORA requer pelo menos ", .MIN_GENES_ORA,
             ". Pulando ORA — considere relaxar os limites.")
  } else {
    ego <- tryCatch(
      enrichGO(
        gene          = genes_sig,
        OrgDb         = org.Hs.eg.db,
        keyType       = "SYMBOL",
        ont           = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff  = padj_ora,
        qvalueCutoff  = padj_ora,
        readable      = TRUE
      ),
      error = function(e) {
        log_error("enrichGO falhou: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(ego) && nrow(ego@result) > 0) {
      # >>>>>>>>>>>>>>>>>>>>>>>
      # Simplificar: remover termos GO altamente redundantes (similaridade semântica > 0.7)
      # <<<<<<<<<<<<<<<<<<<<<<<
      ego <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

      # >>>>>>>>>>>>>>>>>>>>>>>
      # Filtrar: manter apenas termos com >= 2 genes (Count = 1 não é informativo)
      # <<<<<<<<<<<<<<<<<<<<<<<
      ego_filtered <- ego
      ego_filtered@result <- ego@result %>%
        dplyr::filter(Count >= 2)

      n_terms <- nrow(ego_filtered@result)
      log_info("Termos significativos GO ORA (Count ≥ 2, após simplificação): ", n_terms)

      if (n_terms > 0) {
        write_tsv(as.data.frame(ego),
                  file.path(out_dir, paste0("GO_ORA_", target_name, ".tsv")))

        .plot_go_barplot(ego_filtered, target_name, plots_dir)
        .plot_go_dotplot(ego_filtered, target_name, plots_dir)

        results_list$ego <- ego
      } else {
        log_warn("Todos os termos GO ORA têm Count < 2 após filtragem. ",
                 "Resultados salvos em TSV, mas gráficos pulados.")
        write_tsv(as.data.frame(ego),
                  file.path(out_dir, paste0("GO_ORA_", target_name, ".tsv")))
      }
    } else {
      log_warn("Nenhum termo GO significativo encontrado via ORA para: ", target_name)
    }
  }

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 2. Análise de Enriquecimento de Conjunto de Genes (GSEA)
  #    Métrica de ranqueamento: estatística Wald dos resultados brutos do DESeq2
  #    (mais poderoso que log2FC; captura tanto o tamanho do efeito quanto a incerteza)
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("--- Passo 2: GSEA ---")

  gsea_df <- res_raw %>%
    dplyr::filter(!is.na(stat), !is.na(gene_symbol)) %>%
    dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(stat))

  gene_list        <- gsea_df$stat
  names(gene_list) <- gsea_df$gene_symbol

  log_info("Genes na lista ranqueada para GSEA: ", length(gene_list))

  gse <- tryCatch(
    gseGO(
      geneList      = gene_list,
      OrgDb         = org.Hs.eg.db,
      keyType       = "SYMBOL",
      ont           = "BP",
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      pvalueCutoff  = padj_gsea,
      pAdjustMethod = "BH",
      verbose       = FALSE,
      seed          = 42  # para reprodutibilidade
    ),
    error = function(e) {
      log_error("gseGO falhou: ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(gse) && nrow(gse@result) > 0) {
    # >>>>>>>>>>>>>>>>>>>>>>>
    # Simplificar resultados do GSEA
    # <<<<<<<<<<<<<<<<<<<<<<<
    gse <- simplify(gse, cutoff = 0.7, by = "p.adjust", select_fun = min)

    n_terms <- nrow(gse@result)
    log_info("Termos significativos GSEA (após simplificação): ", n_terms)

    write_tsv(as.data.frame(gse),
              file.path(out_dir, paste0("GSEA_", target_name, ".tsv")))

    .plot_gsea_dotplot(gse, target_name, plots_dir)
    .plot_gsea_ridgeplot(gse, target_name, plots_dir)

    results_list$gse <- gse
  } else {
    log_warn("Nenhum termo GSEA significativo encontrado para: ", target_name)
  }

  log_info("Análise de enriquecimento concluída: ", target_name)
  invisible(results_list)
}

# >>>>>>>>>>>>>>>>>>>>>>>
# Auxiliares internos para gráficos
# <<<<<<<<<<<<<<<<<<<<<<<

.wrap_labels <- function(x, width = 40) stringr::str_wrap(x, width = width)

.plot_go_barplot <- function(ego, target_name, plots_dir) {
  p <- barplot(ego, showCategory = 20) +
    labs(
      title    = paste(target_name, "- Barplot de Enriquecimento GO"),
      subtitle = "Processo Biológico | ORA | Simplificado",
      x        = "Contagem de Genes",
      y        = NULL
    ) +
    scale_y_discrete(labels = .wrap_labels) +
    scale_fill_gradient(low = "#8B3A8B", high = "#D4A0D4", name = "FDR") +
    theme_pipeline(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))

  out <- file.path(plots_dir, paste0("GO_Barplot_", target_name, ".png"))
  ggsave(out, plot = p, width = 11, height = 8, dpi = 300)
  log_info("GO Barplot salvo → ", out)
  invisible(p)
}

.plot_go_dotplot <- function(ego, target_name, plots_dir) {
  p <- dotplot(ego, showCategory = 20) +
    labs(
      title    = paste(target_name, "- Dotplot de Enriquecimento GO"),
      subtitle = "Processo Biológico | ORA | Simplificado",
      x        = "Razão de Genes",
      y        = NULL
    ) +
    scale_y_discrete(labels = .wrap_labels) +
    scale_color_gradient(low = "#8B3A8B", high = "#D4A0D4", name = "FDR") +
    theme_pipeline(base_size = 11) +
    theme(axis.text.y = element_text(size = 9))

  out <- file.path(plots_dir, paste0("GO_Dotplot_", target_name, ".png"))
  ggsave(out, plot = p, width = 11, height = 8, dpi = 300)
  log_info("GO Dotplot salvo → ", out)
  invisible(p)
}

.plot_gsea_dotplot <- function(gse, target_name, plots_dir) {
  p <- dotplot(gse, showCategory = 10, split = ".sign") +
    facet_grid(. ~ .sign) +
    labs(
      title    = paste(target_name, "- GSEA Dotplot"),
      subtitle = "Processo Biológico | Ranqueado pela estatística Wald | Simplificado",
      x        = "Razão de Genes",
      y        = NULL
    ) +
    scale_y_discrete(labels = .wrap_labels) +
    scale_color_gradient(
      low    = "firebrick",
      high   = "steelblue",
      name   = "FDR",
      breaks = scales::pretty_breaks(n = 4),
      labels = scales::label_scientific(digits = 1),
      guide  = guide_colorbar(
        barwidth  = 0.8,
        barheight = 5,
        ticks     = TRUE
      )
    ) +
    theme_pipeline(base_size = 11) +
    theme(
      strip.text       = element_text(face = "bold"),
      axis.text.y      = element_text(size = 9),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8),
      legend.key.width = unit(0.4, "cm")
    )

  out <- file.path(plots_dir, paste0("GSEA_Dotplot_", target_name, ".png"))
  ggsave(out, plot = p, width = 14, height = 8, dpi = 300)
  log_info("GSEA Dotplot salvo → ", out)
  invisible(p)
}

.plot_gsea_ridgeplot <- function(gse, target_name, plots_dir) {
  p <- tryCatch(
    ridgeplot(gse, showCategory = 15, fill = "p.adjust") +
      labs(
        title    = paste(target_name, "- GSEA Ridge Plot"),
        subtitle = "Distribuição do ranque dos genes nos conjuntos enriquecidos",
        x        = "Distribuição de Log₂ Fold Change",
        y        = NULL
      ) +
      scale_fill_gradient(low = "firebrick", high = "#FFD700", name = "FDR") +
      scale_y_discrete(labels = .wrap_labels) +
      theme_pipeline(base_size = 11) +
      theme(axis.text.y = element_text(size = 9)),
    error = function(e) {
      log_warn("Ridge plot pulado (poucos termos ou problema no ggridges): ",
               conditionMessage(e))
      NULL
    }
  )

  if (!is.null(p)) {
    out <- file.path(plots_dir, paste0("GSEA_RidgePlot_", target_name, ".png"))
    ggsave(out, plot = p, width = 12, height = 9, dpi = 300)
    log_info("GSEA Ridge plot salvo → ", out)
  }

  invisible(p)
}
