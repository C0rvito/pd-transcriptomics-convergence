// ============================================================================
// Module: FASTP  (adapter trimming + quality filtering)
//
// Output layout:
//   data/reads/trimmed/<GENE>/<SAMPLE>/   trimmed FASTQ.gz
//   data/qc/fastp/<GENE>/<SAMPLE>/        JSON + HTML QC reports
// ============================================================================

process FASTP {
    tag "${meta.gene}/${meta.id}"

    publishDir [
        [
            path: "${params.outdir}/reads/trimmed/${meta.gene}/${meta.id}",
            mode: 'copy',
            pattern: '*.fastq.gz'
        ],
        [
            path: "${params.outdir}/qc/fastp/${meta.gene}/${meta.id}",
            mode: 'copy',
            pattern: '*.{json,html}'
        ]
    ]

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R1.trimmed.fastq.gz"),
                     path("${meta.id}_R2.trimmed.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),           emit: json
    tuple val(meta), path("${meta.id}.fastp.html"),           emit: html

    script:
    """
    fastp \\
        --in1  ${r1} \\
        --in2  ${r2} \\
        --out1 ${meta.id}_R1.trimmed.fastq.gz \\
        --out2 ${meta.id}_R2.trimmed.fastq.gz \\
        --detect_adapter_for_pe \\
        --trim_poly_g \\
        --trim_poly_x \\
        --qualified_quality_phred  ${params.fastp_qual} \\
        --unqualified_percent_limit ${params.fastp_unqual_pct} \\
        --length_required  ${params.fastp_min_len} \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json \\
        --html ${meta.id}.fastp.html
    """
}
