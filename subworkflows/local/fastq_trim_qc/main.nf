//
// FASTQ_TRIM_QC: QC and trimming of raw FASTQ files
//

include { FASTQC } from '../../../modules/nf-core/fastqc/main'
include { FASTP  } from '../../../modules/nf-core/fastp/main'

workflow FASTQ_TRIM_QC {

    take:
    ch_reads  // channel: [ val(meta), [ path(reads) ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: FastQC on raw reads
    //
    FASTQC ( ch_reads )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // MODULE: Fastp for trimming
    //
    ch_reads_for_fastp = ch_reads.map { meta, reads -> [ meta, reads, [] ] }

    FASTP (
        ch_reads_for_fastp,
        false,  // discard_trimmed_pass
        false,  // save_trimmed_fail
        false   // save_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions.first())

    emit:
    reads    = FASTP.out.reads       // channel: [ val(meta), [ path(trimmed_reads) ] ]
    fastqc_html = FASTQC.out.html    // channel: [ val(meta), path(html) ]
    fastqc_zip  = FASTQC.out.zip     // channel: [ val(meta), path(zip) ]
    fastp_json  = FASTP.out.json     // channel: [ val(meta), path(json) ]
    fastp_html  = FASTP.out.html     // channel: [ val(meta), path(html) ]
    fastp_log   = FASTP.out.log      // channel: [ val(meta), path(log) ]
    versions    = ch_versions        // channel: [ path(versions.yml) ]
}
