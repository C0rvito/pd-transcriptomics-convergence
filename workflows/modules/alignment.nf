// ============================================================================
// Module: STAR_ALIGN  (splice-aware alignment + gene counts)
//
// Output layout:
//   data/star/<GENE>/<SAMPLE>/   BAM, ReadsPerGene.out.tab, SJ.out.tab, logs
// ============================================================================

process STAR_ALIGN {
    tag "${meta.gene}/${meta.id}"

    publishDir "${params.outdir}/star/${meta.gene}/${meta.id}",
        mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)
    path  star_index

    output:
    tuple val(meta), path("${meta.id}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.id}.ReadsPerGene.out.tab"),          emit: counts
    tuple val(meta), path("${meta.id}.SJ.out.tab"),                    emit: sj
    tuple val(meta), path("${meta.id}.Log.final.out"),                 emit: log_final
    tuple val(meta), path("${meta.id}.Log.out"),                       emit: log
    tuple val(meta), path("${meta.id}.Log.progress.out"),              emit: log_progress

    script:
    """
    STAR \\
        --runThreadN          ${task.cpus} \\
        --genomeDir           ${star_index} \\
        --readFilesIn         ${r1} ${r2} \\
        --readFilesCommand    gunzip -c \\
        --outFileNamePrefix   ${meta.id}. \\
        --outSAMtype          BAM SortedByCoordinate \\
        --quantMode           GeneCounts \\
        --sjdbOverhang        ${params.star_overhang}
    """
}
