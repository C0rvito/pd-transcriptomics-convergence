#!/usr/bin/env nextflow
// ============================================================================
// Main workflow – RNA-seq preprocessing and alignment
//
// Input layout:
//   params.reads_dir/<GENE>/files/<SRR>_{1,2}.fastq.gz
//
// Output layout (all under params.outdir = data/):
//   qc/raw/<GENE>/<SAMPLE>/          FastQC on raw reads
//   qc/trimmed/<GENE>/<SAMPLE>/      FastQC on trimmed reads
//   qc/fastp/<GENE>/<SAMPLE>/        fastp JSON + HTML reports
//   qc/multiqc/<GENE>/               MultiQC aggregated report per gene
//   reads/trimmed/<GENE>/<SAMPLE>/   Trimmed FASTQ.gz files
//   star/<GENE>/<SAMPLE>/            STAR BAM + counts + logs
// ============================================================================

nextflow.enable.dsl = 2

include { FASTQC as FASTQC_RAW     } from './modules/fastqc'
include { FASTQC as FASTQC_TRIMMED } from './modules/fastqc'
include { FASTP                    } from './modules/trimming'
include { STAR_ALIGN               } from './modules/alignment'
include { MULTIQC                  } from './modules/multiqc'

// ----------------------------------------------------------------------------
workflow {

    // Build input channel: [meta, r1, r2]
    // meta = [id: <SRR>, gene: <GENE>]
    reads_ch = Channel
        .fromFilePairs(
            "${params.reads_dir}/*/files/*_{1,2}.fastq.gz",
            checkIfExists: true
        )
        .map { sample_id, files ->
            // files[0] path: .../raw/<GENE>/files/<SRR>_1.fastq.gz
            def gene = files[0].parent.parent.name
            def meta = [id: sample_id, gene: gene]
            tuple(meta, files[0], files[1])
        }

    // 1 ── FastQC on raw reads ────────────────────────────────────────────────
    FASTQC_RAW(
        reads_ch.map { meta, r1, r2 -> [meta + [step: 'raw'], r1, r2] }
    )

    // 2 ── Trimming with fastp ────────────────────────────────────────────────
    FASTP(reads_ch)

    // 3 ── FastQC on trimmed reads ────────────────────────────────────────────
    FASTQC_TRIMMED(
        FASTP.out.reads.map { meta, r1, r2 -> [meta + [step: 'trimmed'], r1, r2] }
    )

    // 4 ── STAR alignment ─────────────────────────────────────────────────────
    star_index_ch = Channel.value(file(params.star_index, checkIfExists: true))
    STAR_ALIGN(FASTP.out.reads, star_index_ch)

    // 5 ── MultiQC per gene ───────────────────────────────────────────────────
    // Collect every QC artefact, tag by gene, then group for a per-gene report
    qc_files_ch = Channel.empty()
        .mix(
            FASTQC_RAW.out.zip.map       { meta, f -> [meta.gene, f] },
            FASTQC_RAW.out.html.map      { meta, f -> [meta.gene, f] },
            FASTP.out.json.map           { meta, f -> [meta.gene, f] },
            FASTP.out.html.map           { meta, f -> [meta.gene, f] },
            FASTQC_TRIMMED.out.zip.map   { meta, f -> [meta.gene, f] },
            STAR_ALIGN.out.log_final.map { meta, f -> [meta.gene, f] }
        )
        .groupTuple()
        .map { gene, files -> tuple(gene, files.flatten()) }

    MULTIQC(qc_files_ch)
}
