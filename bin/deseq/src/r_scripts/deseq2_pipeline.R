# >>>>>>>>>>>>>>>>>>>>>>>
# Arquivo: src/deseq2_pipeline.R
# Descrição: Análise de Expressão Diferencial DESeq2
#   - Encolhimento de LFC via apeglm
#   - Retorna resultados encolhidos (volcano/MA) e estatística bruta (ranking GSEA)
#   - Gera gráficos de CQ: PCA, dispersão, fatores de tamanho
# <<<<<<<<<<<<<<<<<<<<<<<

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(readr)
  library(tibble)
  library(apeglm)
  library(ggplot2)
})

source("src/utils.R")

# >>>>>>>>>>>>>>>>>>>>>>>
#' Executa a análise de expressão diferencial DESeq2
#'
#' @param target_name  Texto. Rótulo do coorte/dataset (ex: "LRRK2").
#' @param counts_file  Caminho para o CSV com colunas: gene_id, gene_symbol, <amostras>.
#' @param meta_file    Caminho para o CSV com colunas: sample, condition, [covariáveis].
#' @param lfc_threshold |log2FC| mínimo para definir um DEG significativo (padrão 1).
#' @param padj_threshold Limite de FDR (padrão 0.05).
#' @param min_count    Contagem mínima de leituras para pré-filtragem (padrão 10).
#' @param out_dir      Diretório de saída para os resultados (padrão "results").
#' @param plots_dir    Diretório de saída para gráficos de CQ (padrão "plots").
#'
#' @return Lista nomeada: dds, res_shrunken, res_raw, norm_counts, degs
# <<<<<<<<<<<<<<<<<<<<<<<
run_deseq2_analysis <- function(
    target_name,
    counts_file,
    meta_file,
    lfc_threshold  = 1,
    padj_threshold = 0.05,
    min_count      = 10,
    out_dir        = "results",
    plots_dir      = "plots"
) {
  log_info("=== Análise DESeq2: ", target_name, " ===")
  setup_dirs(c(out_dir, plots_dir))
  assert_file(counts_file, "Arquivo de contagens")
  assert_file(meta_file,   "Arquivo de metadados")

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 1. Carregar e preparar metadados
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("Carregando metadados...")
  metadata <- read_csv(meta_file, show_col_types = FALSE)

  required_cols <- c("sample", "condition")
  missing_cols  <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Metadados faltando colunas obrigatórias: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  meta_clean <- metadata %>%
    dplyr::mutate(condition = factor(condition, levels = c("Control", "Mutant"))) %>%
    as.data.frame()
  rownames(meta_clean) <- meta_clean$sample

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 2. Carregar e preparar matriz de contagem
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("Carregando matriz de contagens...")
  counts_df <- read_csv(counts_file, show_col_types = FALSE)

  required_count_cols <- c("gene_id", "gene_symbol")
  missing_count_cols  <- setdiff(required_count_cols, colnames(counts_df))
  if (length(missing_count_cols) > 0) {
    stop("Arquivo de contagens faltando colunas obrigatórias: ",
         paste(missing_count_cols, collapse = ", "), call. = FALSE)
  }

  gene_annotations <- counts_df %>%
    dplyr::select(gene_id, gene_symbol) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE)

  counts_mat <- counts_df %>%
    dplyr::select(-gene_symbol) %>%
    tibble::column_to_rownames("gene_id") %>%
    as.matrix()

  # >>>>>>>>>>>>>>>>>>>>>>>
  # Garantir que as contagens sejam inteiras
  # <<<<<<<<<<<<<<<<<<<<<<<
  storage.mode(counts_mat) <- "integer"

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 3. Alinhar amostras entre contagens e metadados
  # <<<<<<<<<<<<<<<<<<<<<<<
  shared_samples <- intersect(colnames(counts_mat), rownames(meta_clean))

  if (length(shared_samples) < 3) {
    stop("Menos de 3 amostras compartilhadas entre contagens e metadados. ",
         "Verifique a consistência dos nomes das amostras.", call. = FALSE)
  }

  if (length(shared_samples) < ncol(counts_mat)) {
    dropped <- setdiff(colnames(counts_mat), shared_samples)
    log_warn(length(dropped), " coluna(s) de contagem não encontradas nos metadados e serão removidas: ",
             paste(dropped, collapse = ", "))
  }

  meta_aligned   <- meta_clean[shared_samples, , drop = FALSE]
  counts_aligned <- counts_mat[, shared_samples, drop = FALSE]
  stopifnot(identical(colnames(counts_aligned), rownames(meta_aligned)))

  log_info("Amostras por condição:")
  print(table(meta_aligned$condition))

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 4. Construir DESeqDataSet e pré-filtrar genes de baixa contagem
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("Construindo DESeqDataSet...")
  dds <- DESeqDataSetFromMatrix(
    countData = counts_aligned,
    colData   = meta_aligned,
    design    = ~ condition
  )

  n_before <- nrow(dds)
  min_samples <- min(table(meta_aligned$condition))
  keep        <- rowSums(counts(dds) >= min_count) >= min_samples
  dds         <- dds[keep, ]
  log_info("Pré-filtragem: ", n_before, " → ", nrow(dds),
           " genes (mínimo ", min_count, " contagens em ≥", min_samples, " amostras)")

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 5. Executar DESeq2
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("Ajustando modelo DESeq2...")
  dds <- DESeq(dds, quiet = TRUE)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # CQ: Gráfico de dispersão
  # <<<<<<<<<<<<<<<<<<<<<<<
  .save_dispersion_plot(dds, target_name, plots_dir)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # CQ: PCA em contagens transformadas por rlog
  # <<<<<<<<<<<<<<<<<<<<<<<
  .save_pca_plot(dds, meta_aligned, target_name, plots_dir)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # CQ: Fatores de tamanho
  # <<<<<<<<<<<<<<<<<<<<<<<
  .save_size_factor_plot(dds, meta_aligned, target_name, plots_dir)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 6. Extrair resultados
  #    - res_raw:      Estatística Wald disponível → usada para ranking GSEA
  #    - res_shrunken: Encolhimento LFC apeglm → usado para gráficos volcano / MA
  # <<<<<<<<<<<<<<<<<<<<<<<
  coef_name <- "condition_Mutant_vs_Control"
  log_info("Extraindo resultados brutos (estatística Wald)...")
  res_raw <- results(dds, name = coef_name, independentFiltering = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotations, by = "gene_id") %>%
    dplyr::relocate(gene_symbol, .after = gene_id) %>%
    dplyr::arrange(padj)

  log_info("Aplicando encolhimento de LFC (apeglm)...")
  res_shrunken <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotations, by = "gene_id") %>%
    dplyr::relocate(gene_symbol, .after = gene_id) %>%
    dplyr::arrange(padj)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 7. Definir DEGs significativos
  #    padj < threshold E |log2FC| >= lfc_threshold (das estimativas encolhidas)
  # <<<<<<<<<<<<<<<<<<<<<<<
  degs <- res_shrunken %>%
    dplyr::filter(padj < padj_threshold, abs(log2FoldChange) >= lfc_threshold)

  n_up   <- sum(degs$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(degs$log2FoldChange < 0, na.rm = TRUE)
  log_info("DEGs significativos: ", nrow(degs),
           " (aumentados: ", n_up, ", reduzidos: ", n_down, ")")

  if (nrow(degs) < 5) {
    log_warn("Poucos DEGs detectados (n = ", nrow(degs), "). ",
             "Considere relaxar lfc_threshold ou padj_threshold para análise de enriquecimento.")
  }

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 8. Contagens normalizadas
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("Extraindo contagens normalizadas...")
  norm_counts <- counts(dds, normalized = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::left_join(gene_annotations, by = "gene_id") %>%
    dplyr::relocate(gene_symbol, .after = gene_id)

  # >>>>>>>>>>>>>>>>>>>>>>>
  # 9. Exportar todas as tabelas de resultados
  # <<<<<<<<<<<<<<<<<<<<<<<
  log_info("Gravando tabelas de resultados em: ", out_dir)
  write_tsv(res_shrunken,
            file.path(out_dir, paste0("DE_Mutant_vs_Control_", target_name, ".tsv")))
  write_tsv(res_raw,
            file.path(out_dir, paste0("DE_Raw_Stat_",          target_name, ".tsv")))
  write_tsv(norm_counts,
            file.path(out_dir, paste0("NormCounts_",            target_name, ".tsv")))
  write_tsv(degs,
            file.path(out_dir, paste0("DEGs_",                  target_name, ".tsv")))

  log_info("Análise DESeq2 concluída: ", target_name)

  invisible(list(
    dds          = dds,
    res_shrunken = res_shrunken,
    res_raw      = res_raw,
    norm_counts  = norm_counts,
    degs         = degs
  ))
}

# >>>>>>>>>>>>>>>>>>>>>>>
# Auxiliares internos para gráficos de CQ (prefixados com . para indicar privado)
# <<<<<<<<<<<<<<<<<<<<<<<

.save_dispersion_plot <- function(dds, target_name, plots_dir) {
  out <- file.path(plots_dir, paste0("QC_Dispersion_", target_name, ".png"))
  png(out, width = 800, height = 600, res = 150)
  DESeq2::plotDispEsts(dds,
    main = paste("Estimativas de Dispersão —", target_name))
  dev.off()
  log_info("CQ: Gráfico de dispersão salvo → ", out)
}

.save_pca_plot <- function(dds, meta_aligned, target_name, plots_dir) {
  rld    <- rlog(dds, blind = TRUE)
  pca_data <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

  p <- ggplot(pca_data, aes(PC1, PC2, color = condition, label = name)) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Control" = "steelblue", "Mutant" = "firebrick")) +
    labs(
      title    = paste("PCA — contagens rlog:", target_name),
      subtitle = "Mutant vs Control iPSCs",
      x = paste0("PC1: ", pct_var[1], "% de variância"),
      y = paste0("PC2: ", pct_var[2], "% de variância")
    ) +
    theme_pipeline()

  out <- file.path(plots_dir, paste0("QC_PCA_", target_name, ".png"))
  ggsave(out, plot = p, width = 7, height = 5, dpi = 300)
  log_info("CQ: Gráfico PCA salvo → ", out)
}

.save_size_factor_plot <- function(dds, meta_aligned, target_name, plots_dir) {
  sf_df <- data.frame(
    sample    = colnames(dds),
    size_factor = sizeFactors(dds),
    condition  = meta_aligned$condition
  )

  p <- ggplot(sf_df, aes(x = reorder(sample, size_factor),
                          y = size_factor, fill = condition)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("Control" = "steelblue", "Mutant" = "firebrick")) +
    labs(
      title = paste("Fatores de Tamanho —", target_name),
      x = "Amostra", y = "Fator de Tamanho"
    ) +
    theme_pipeline() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  out <- file.path(plots_dir, paste0("QC_SizeFactors_", target_name, ".png"))
  ggsave(out, plot = p, width = 8, height = 4, dpi = 300)
  log_info("CQ: Gráfico de fatores de tamanho salvo → ", out)
}
