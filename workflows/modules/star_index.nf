// ============================================================================
// Module: STAR_INDEX  (genome index generation – run once, use many times)
//
// Usage in main.nf (optional step):
//
//   include { STAR_INDEX } from './modules/star_index'
//
//   STAR_INDEX(
//       file(params.genome_fasta),
//       file(params.gtf),
//       params.star_overhang
//   )
//   // Pass STAR_INDEX.out.index as star_index_ch to STAR_ALIGN
//
// Output:
//   params.star_index/   (written directly to the configured index directory)
// ============================================================================

process STAR_INDEX {
    tag "genome index"

    publishDir "${params.star_index}",
        mode: 'copy',
        overwrite: true

    input:
    path  genome_fasta
    path  gtf
    val   overhang

    output:
    path "STAR_Index/", emit: index

    script:
    """
    mkdir -p STAR_Index

    STAR \\
        --runThreadN      ${task.cpus} \\
        --runMode         genomeGenerate \\
        --genomeDir       STAR_Index \\
        --genomeFastaFiles ${genome_fasta} \\
        --sjdbGTFfile     ${gtf} \\
        --sjdbOverhang    ${overhang}
    """
}
