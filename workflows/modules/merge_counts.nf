// ============================================================================
// Module: MERGE_COUNTS
//
// Merge STAR ReadsPerGene.out.tab files into one count matrix per gene.
// ============================================================================

process MERGE_COUNTS {
    label 'process_low'

    publishDir { "${params.outdir}/de/counts/${gene}" }, mode: 'copy'

    input:
    tuple val(gene), path(gene_dir)

    output:
    path "star_counts_brutas.tsv", emit: counts

    script:
    """
    python3 ${params.project_root}/bin/merge_star_counts.py \\
        --star-dir . \\
        --outdir .
    
    # The script creates <outdir>/<GENE>/star_counts_brutas.tsv
    # We move it to the current directory for Nextflow to pick it up easily
    mv ${gene}/star_counts_brutas.tsv .
    """
}
