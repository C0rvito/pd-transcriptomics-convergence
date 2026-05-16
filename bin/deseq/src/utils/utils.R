# >>>>>>>>>>>>>>>>>>>>>>>
# Arquivo: src/utils.R
# Descrição: Funções utilitárias — logs, configuração de diretórios, exportação de sessão
# <<<<<<<<<<<<<<<<<<<<<<<

# >>>>>>>>>>>>>>>>>>>>>>>
#' Log de INFORMAÇÃO com carimbo de tempo
# <<<<<<<<<<<<<<<<<<<<<<<
log_info <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S]"), "INFO  |", paste0(...), "\n")
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Log de AVISO com carimbo de tempo
# <<<<<<<<<<<<<<<<<<<<<<<
log_warn <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S]"), "WARN  |", paste0(...), "\n")
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Log de ERRO com carimbo de tempo (não interrompe a execução)
# <<<<<<<<<<<<<<<<<<<<<<<
log_error <- function(...) {
  cat(format(Sys.time(), "[%H:%M:%S]"), "ERROR |", paste0(...), "\n")
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Criar vários diretórios de saída de uma vez
# <<<<<<<<<<<<<<<<<<<<<<<
setup_dirs <- function(dirs) {
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Verificar se um arquivo existe ou lançar um erro descritivo
# <<<<<<<<<<<<<<<<<<<<<<<
assert_file <- function(path, label = "File") {
  if (!file.exists(path)) {
    stop(label, " não encontrado: ", path, call. = FALSE)
  }
  invisible(path)
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Salvar informações da sessão R em results/ para reprodutibilidade
# <<<<<<<<<<<<<<<<<<<<<<<
save_session_info <- function(out_dir = "results") {
  out_path <- file.path(out_dir, "session_info.txt")
  writeLines(capture.output(utils::sessionInfo()), out_path)
  log_info("Informações da sessão salvas em: ", out_path)
}

# >>>>>>>>>>>>>>>>>>>>>>>
#' Tema ggplot2 consistente para o pipeline
# <<<<<<<<<<<<<<<<<<<<<<<
theme_pipeline <- function(base_size = 13) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(color = "grey40", size = base_size - 1),
      legend.position = "bottom",
      axis.text     = ggplot2::element_text(color = "black")
    )
}
