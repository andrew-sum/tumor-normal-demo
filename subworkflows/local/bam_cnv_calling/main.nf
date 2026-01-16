/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_CNV_CALLING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Copy number variation calling with CNVkit
    - Pairs tumor samples with matched normal by patient ID
    - Runs CNVkit batch in tumor-normal mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNVKIT_BATCH } from '../../../modules/nf-core/cnvkit/batch/main'

workflow BAM_CNV_CALLING {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fai              // channel: [ val(meta), path(fai) ]
    ch_targets          // channel: [ path(bed) ] - optional

    main:
    ch_versions = Channel.empty()

    //
    // Separate tumor and normal samples by status
    // status: 0 = normal, 1 = tumor
    //
    ch_bam_bai
        .branch { meta, bam, bai ->
            tumor:  meta.status == 1
                return [ meta, bam, bai ]
            normal: meta.status == 0
                return [ meta, bam, bai ]
        }
        .set { ch_branched }

    //
    // Pair tumor with matched normal by patient ID
    // Output: [ meta_tumor, tumor_bam, normal_bam ]
    //
    ch_branched.tumor
        .map { meta, bam, bai -> [ meta.patient, meta, bam, bai ] }
        .set { ch_tumor_by_patient }

    ch_branched.normal
        .map { meta, bam, bai -> [ meta.patient, meta, bam, bai ] }
        .set { ch_normal_by_patient }

    ch_tumor_by_patient
        .join(ch_normal_by_patient)
        .map { patient, tumor_meta, tumor_bam, tumor_bai, normal_meta, normal_bam, normal_bai ->
            // Create new meta with both sample names
            def new_meta = tumor_meta.clone()
            new_meta.normal_id = normal_meta.id
            [ new_meta, tumor_bam, normal_bam ]
        }
        .set { ch_paired }

    //
    // Prepare target regions channel - wrap in tuple with empty meta
    //
    ch_targets_input = ch_targets
        .map { targets -> [ [:], targets ] }
        .ifEmpty { [ [:], [] ] }

    //
    // MODULE: CNVKIT_BATCH
    // Input: tuple val(meta), path(tumor), path(normal)
    //        tuple val(meta2), path(fasta)
    //        tuple val(meta3), path(fasta_fai)
    //        tuple val(meta4), path(targets)
    //        tuple val(meta5), path(reference)
    //        val panel_of_normals
    //
    CNVKIT_BATCH (
        ch_paired,
        ch_fasta,
        ch_fai,
        ch_targets_input,
        Channel.value([ [:], [] ]),    // reference (empty - will be computed)
        false                          // panel_of_normals
    )
    ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions.first())

    emit:
    cns      = CNVKIT_BATCH.out.cns      // channel: [ val(meta), path(cns) ]
    cnr      = CNVKIT_BATCH.out.cnr      // channel: [ val(meta), path(cnr) ]
    pdf      = CNVKIT_BATCH.out.pdf      // channel: [ val(meta), path(pdf) ]
    bed      = CNVKIT_BATCH.out.bed      // channel: [ val(meta), path(bed) ]
    versions = ch_versions               // channel: [ path(versions.yml) ]
}
