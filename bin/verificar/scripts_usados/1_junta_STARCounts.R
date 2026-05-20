# 1. Definir o diretório onde estão os arquivos .tab
directory <- "/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/star_teste"
files <- list.files(path = directory, pattern = "\\.ReadsPerGene\\.out\\.tab$", full.names = TRUE)

# 2. Função para ler apenas a 1ª coluna (ID) e a 4ª coluna (Counts)
read_counts <- function(file) {
  # Lendo o arquivo pulando as 4 primeiras linhas (metadados do STAR)
  data <- read.table(file, skip = 4, header = FALSE)
  # Coluna 1 = GeneID, Coluna 4 = Reverse Counts
  return(data[, c(1, 4)])
}

# 3. Carregar a primeira amostra para iniciar a matriz
first_sample_name <- gsub(".ReadsPerGene.out.tab", "", basename(files[1]))
master_matrix <- read_counts(files[1])
colnames(master_matrix) <- c("GeneID", first_sample_name)

# 4. Loop para juntar as demais amostras
for (i in 2:length(files)) {
  sample_name <- gsub(".ReadsPerGene.out.tab", "", basename(files[i]))
  current_sample <- read_counts(files[i])
  colnames(current_sample) <- c("GeneID", sample_name)
  
  # Faz o merge baseado no ID do Gene
  master_matrix <- merge(master_matrix, current_sample, by = "GeneID")
}

# 5. Definir o GeneID como nome das linhas e remover a coluna repetida
rownames(master_matrix) <- master_matrix$GeneID
master_matrix$GeneID <- NULL

# 6. Salvar a matriz final para usar no DESeq2
write.csv(master_matrix, "/home/jess/Documents/estudos/pd-transcriptomics-convergence/data/GSE149632/matriz_contagem_final.csv")

print("Matriz consolidada com sucesso!")
