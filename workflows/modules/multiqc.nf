// ============================================================================
// Module: MULTIQC  (per-gene aggregated QC report)
//
// Receives all QC artefacts for a single gene (FastQC zips/HTMLs, fastp
// JSONs, STAR Log.final.out files) and generates one MultiQC report.
//
// Output layout:
//   data/qc/multiqc/<GENE>/
// ============================================================================

process MULTIQC {
    tag "${gene}"

    publishDir "${params.outdir}/qc/multiqc/${gene}",
        mode: 'copy'

    input:
    tuple val(gene), path(qc_files)

    output:
    tuple val(gene), path("*.html"),  emit: report
    tuple val(gene), path("*_data/"), emit: data, optional: true

    script:
    """
    multiqc --title "QC Report ${gene}" . -o .
    """
}
