// ============================================================================
// Module: FASTQC
//
// Runs FastQC on a paired-end sample (R1 + R2 in a single call).
// meta.step must be set to 'raw' or 'trimmed' by the caller so that outputs
// land in the correct sub-directory:
//   data/qc/raw/<GENE>/<SAMPLE>/
//   data/qc/trimmed/<GENE>/<SAMPLE>/
// ============================================================================

process FASTQC {
    tag "${meta.gene}/${meta.id} [${meta.step}]"

    publishDir "${params.outdir}/qc/${meta.step}/${meta.gene}/${meta.id}",
        mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*.zip"),  emit: zip
    tuple val(meta), path("*.html"), emit: html

    script:
    """
    fastqc -t ${task.cpus} -o . ${r1} ${r2}
    """
}
