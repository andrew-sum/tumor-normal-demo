//
// FASTQ_ALIGN: Align reads to reference genome using BWA-MEM2
//

include { BWAMEM2_MEM    } from '../../../modules/nf-core/bwamem2/mem/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN {

    take:
    ch_reads    // channel: [ val(meta), [ path(reads) ] ]
    ch_fasta    // channel: [ val(meta), path(fasta) ]
    ch_bwamem2  // channel: [ val(meta), path(index) ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: BWA-MEM2 alignment
    // Using sort_bam=true to output coordinate-sorted BAM
    //
    BWAMEM2_MEM (
        ch_reads,
        ch_bwamem2,
        ch_fasta,
        true  // sort_bam - outputs sorted BAM directly
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions.first())

    //
    // MODULE: Index BAM files
    //
    SAMTOOLS_INDEX ( BWAMEM2_MEM.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    // Join BAM and BAI for downstream tools
    ch_bam_bai = BWAMEM2_MEM.out.bam.join(SAMTOOLS_INDEX.out.bai)

    emit:
    bam      = BWAMEM2_MEM.out.bam     // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai  // channel: [ val(meta), path(bai) ]
    bam_bai  = ch_bam_bai              // channel: [ val(meta), path(bam), path(bai) ]
    versions = ch_versions             // channel: [ path(versions.yml) ]
}
