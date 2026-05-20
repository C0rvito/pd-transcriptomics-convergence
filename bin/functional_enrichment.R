#!/usr/bin/env Rscript

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(openxlsx)

# Usage: Rscript functional_enrichment.R <de_results.csv> <outdir> <comparison_name>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript functional_enrichment.R <de_results.csv> <outdir> <comparison_name>")
}

de_file <- args[1]
outdir <- args[2]
comp_name <- args[3]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 1. Load DE Results
de_res <- read.csv(de_file, row.names = 1)

# Filter significant genes (padj < 0.05)
sig_genes_all <- de_res[which(de_res$padj < 0.05), ]
sig_genes_up <- sig_genes_all[which(sig_genes_all$log2FoldChange > 0), ]
sig_genes_down <- sig_genes_all[which(sig_genes_all$log2FoldChange < 0), ]

run_enrichment <- function(genes_df, prefix) {
    if (nrow(genes_df) < 10) {
        message(paste("Too few genes for", prefix, "skipping."))
        return(NULL)
    }
    
    # Remove version from Ensembl IDs if present (e.g., ENSG000001.1 -> ENSG000001)
    gene_ids <- gsub("\\..*", "", rownames(genes_df))
    
    # ── GO Enrichment ────────────────────────────────────────────────────────
    ego <- enrichGO(gene          = gene_ids,
                    OrgDb         = org.Hs.eg.db,
                    keyType       = 'ENSEMBL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.05,
                    qvalueCutoff  = 0.2,
                    readable      = TRUE)
    
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        # Dotplot
        p <- dotplot(ego, showCategory=20, split="ONTOLOGY") + 
             facet_grid(ONTOLOGY~., scale="free") +
             labs(title = paste("GO Enrichment -", comp_name, "(", prefix, ")"))
        ggsave(file.path(outdir, paste0("GO_dotplot_", prefix, ".png")), p, width=10, height=12, dpi=300)
        
        # Save table
        write.xlsx(as.data.frame(ego), file.path(outdir, paste0("GO_results_", prefix, ".xlsx")))
    }
    
    # ── KEGG Enrichment ──────────────────────────────────────────────────────
    # KEGG needs Entrez IDs
    entrez_ids <- mapIds(org.Hs.eg.db, keys=gene_ids, column="ENTREZID", keytype="ENSEMBL")
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    kk <- enrichKEGG(gene         = entrez_ids,
                     organism     = 'hsa',
                     pvalueCutoff = 0.05)
    
    if (!is.null(kk) && nrow(as.data.frame(kk)) > 0) {
        p <- dotplot(kk, showCategory=20) + 
             labs(title = paste("KEGG Enrichment -", comp_name, "(", prefix, ")"))
        ggsave(file.path(outdir, paste0("KEGG_dotplot_", prefix, ".png")), p, width=10, height=8, dpi=300)
        
        # Save table
        write.xlsx(as.data.frame(kk), file.path(outdir, paste0("KEGG_results_", prefix, ".xlsx")))
    }
}

# Run for Up, Down and All
message("Running enrichment for All significant genes...")
run_enrichment(sig_genes_all, "All")
message("Running enrichment for Up-regulated genes...")
run_enrichment(sig_genes_up, "Up")
message("Running enrichment for Down-regulated genes...")
run_enrichment(sig_genes_down, "Down")

message(paste("Enrichment complete for", comp_name))
