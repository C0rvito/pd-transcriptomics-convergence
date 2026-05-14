# >>>>>>>>>>>>>>>>>>>>>>>
# Arquivo: main.R
# Descrição: Orquestrador do pipeline para DE de RNA-Seq + Análise de Enriquecimento
#
# Uso:
#   Rscript main.R
#   OU source("main.R") de uma sessão R
#
# Arquivos de entrada esperados:
#   data/GSE90469/counts/LRRK2_mRNA_counts_with_symbols.csv
#     Colunas: gene_id, gene_symbol, <amostra_1>, <amostra_2>, ...
#   data/GSE90469/counts/metadata.csv
#     Colunas: sample, condition (Control | Mutant)
#
# Estrutura de saída:
#   results/   ← Tabelas TSV (resultados DE, DEGs, contagens normalizadas, enriquecimento)
#   plots/     ← Figuras PNG (CQ, volcano, MA, GO, GSEA)
# <<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>
# 0. Inicialização
# <<<<<<<<<<<<<<<<<<<<<<<
cat(rep("=", 60), "\n", sep = "")
cat("  Pipeline de DE de RNA-Seq + Enriquecimento\n")
cat(format(Sys.time(), "  Iniciado em: %Y-%m-%d %H:%M:%S\n"))
cat(rep("=", 60), "\n", sep = "")

source("src/utils.R")
source("src/deseq2_pipeline.R")
source("src/visualization.R")
source("src/gsea_pipeline.R")

# >>>>>>>>>>>>>>>>>>>>>>>
# Resolver conflitos de namespace introduzidos por pacotes do Bioconductor
# (AnnotationDbi::select e stats::filter mascaram equivalentes do dplyr)
# <<<<<<<<<<<<<<<<<<<<<<<
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

setup_dirs(c("results", "plots"))

# >>>>>>>>>>>>>>>>>>>>>>>
# 1. Configuração da análise
# <<<<<<<<<<<<<<<<<<<<<<<
# Para adicionar mais coortes, anexe entradas a esta lista.
# Cada entrada aciona uma execução completa de DE + enriquecimento.
cohorts <- list(
  LRRK2 = list(
    counts_file = "data/GSE90469/counts/LRRK2_mRNA_counts_with_symbols.csv",
    meta_file   = "data/GSE90469/counts/metadata.csv"
  )
)

# >>>>>>>>>>>>>>>>>>>>>>>
# Limites globais — aplicados consistentemente em todos os coortes
# <<<<<<<<<<<<<<<<<<<<<<<
thresholds <- list(
  lfc_threshold  = 1,      # |log2FC| ≥ 1 para chamar um DEG
  padj_threshold = 0.05,   # FDR < 0.05
  min_count      = 10,     # pré-filtragem: contagem bruta mínima por gene
  padj_ora       = 0.05,   # corte de FDR para GO ORA
  padj_gsea      = 0.05,   # corte de FDR para GSEA
  min_gs_size    = 15,     # tamanho mínimo do conjunto de genes para GSEA
  max_gs_size    = 500     # tamanho máximo do conjunto de genes para GSEA
)

# >>>>>>>>>>>>>>>>>>>>>>>
# 2. Executar pipeline para cada coorte
# <<<<<<<<<<<<<<<<<<<<<<<
pipeline_status <- lapply(names(cohorts), function(target_name) {
  cfg <- cohorts[[target_name]]

  cat("\n", rep("-", 60), "\n", sep = "")
  log_info("Processando coorte: ", target_name)
  cat(rep("-", 60), "\n", sep = "")

  tryCatch({
    # Passo 1 — Expressão Diferencial (DESeq2)
    de_results <- run_deseq2_analysis(
      target_name    = target_name,
      counts_file    = cfg$counts_file,
      meta_file      = cfg$meta_file,
      lfc_threshold  = thresholds$lfc_threshold,
      padj_threshold = thresholds$padj_threshold,
      min_count      = thresholds$min_count
    )

    # Passo 2 — Visualização (Volcano + MA)
    generate_volcano_plot(
      res_shrunken   = de_results$res_shrunken,
      target_name    = target_name,
      lfc_threshold  = thresholds$lfc_threshold,
      padj_threshold = thresholds$padj_threshold
    )

    generate_ma_plot(
      res_shrunken   = de_results$res_shrunken,
      target_name    = target_name,
      padj_threshold = thresholds$padj_threshold
    )

    # Passo 3 — Enriquecimento (GO ORA + GSEA)
    run_enrichment_analysis(
      res_shrunken = de_results$res_shrunken,
      res_raw      = de_results$res_raw,
      degs         = de_results$degs,
      target_name  = target_name,
      padj_ora     = thresholds$padj_ora,
      padj_gsea    = thresholds$padj_gsea,
      min_gs_size  = thresholds$min_gs_size,
      max_gs_size  = thresholds$max_gs_size
    )

    log_info("Coorte concluído: ", target_name)
    list(cohort = target_name, status = "SUCCESS", error = NA)

  }, error = function(e) {
    log_error("Coorte FALHOU: ", target_name)
    log_error("Razão: ", conditionMessage(e))
    list(cohort = target_name, status = "FAILED", error = conditionMessage(e))
  })
})

# >>>>>>>>>>>>>>>>>>>>>>>
# 3. Relatório de resumo
# <<<<<<<<<<<<<<<<<<<<<<<
cat("\n", rep("=", 60), "\n", sep = "")
cat("  Resumo do Pipeline\n")
cat(rep("=", 60), "\n", sep = "")

for (s in pipeline_status) {
  status_icon <- if (s$status == "SUCCESS") "✓" else "✗"
  cat(sprintf("  %s  %-20s %s\n", status_icon, s$cohort, s$status))
  if (!is.na(s$error)) cat("       Erro:", s$error, "\n")
}

save_session_info(out_dir = "results")
cat(format(Sys.time(), "\n  Finalizado em: %Y-%m-%d %H:%M:%S\n"))
cat(rep("=", 60), "\n", sep = "")
