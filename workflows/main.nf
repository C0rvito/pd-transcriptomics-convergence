#!/usr/bin/env nextflow
// ============================================================================
// Main workflow – RNA-seq preprocessing and alignment
// ============================================================================

nextflow.enable.dsl = 2

// include { FASTQC as FASTQC_RAW     } from './modules/fastqc'
// include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
// include { FASTP                    } from './modules/trimming'
// include { STAR_ALIGN               } from './modules/alignment'
// include { MULTIQC                  } from './modules/multiqc'
include { MERGE_COUNTS             } from './modules/merge_counts'

// ----------------------------------------------------------------------------
workflow {

    // // Build input channel: [meta, r1, r2]
    // reads_ch = Channel
    //     .fromFilePairs(
    //         "${params.reads_dir}/*/files/*_{1,2}.fastq.gz",
    //         checkIfExists: true
    //     )
    //     .map { sample_id, files ->
    //         def gene = files[0].parent.parent.name
    //         def meta = [id: sample_id, gene: gene]
    //         tuple(meta, files[0], files[1])
    //     }

    // FASTQC_RAW(
    //     reads_ch.map { meta, r1, r2 -> [meta + [step: 'raw'], r1, r2] }
    // )
    // FASTP(reads_ch)
    // FASTQC_TRIMMED(
    //     FASTP.out.reads.map { meta, r1, r2 -> [meta + [step: 'trimmed'], r1, r2] }
    // )
    // star_index_ch = Channel.value(file(params.star_index, checkIfExists: true))
    // STAR_ALIGN(FASTP.out.reads, star_index_ch)

    // qc_files_ch = Channel.empty()
    //     .mix(
    //         FASTQC_RAW.out.zip.map       { meta, f -> [meta.gene, f] },
    //         FASTQC_RAW.out.html.map      { meta, f -> [meta.gene, f] },
    //         FASTP.out.json.map           { meta, f -> [meta.gene, f] },
    //         FASTP.out.html.map           { meta, f -> [meta.gene, f] },
    //         FASTQC_TRIMMED.out.zip.map   { meta, f -> [meta.gene, f] },
    //         STAR_ALIGN.out.log_final.map { meta, f -> [meta.gene, f] }
    //     )
    //     .groupTuple()
    //     .map { gene, files -> tuple(gene, files.flatten()) }

    // MULTIQC(qc_files_ch)

    // ── MERGE COUNTS ONLY ────────────────────────────────────────────────────
    // Create a channel from existing STAR output directories
    star_dirs_ch = Channel
        .fromPath("${params.outdir}/star/*", type: 'dir')
        .filter { it.name != 'LRRK2' }
        .map { dir -> tuple(dir.name, dir) }

    MERGE_COUNTS(star_dirs_ch)
}
