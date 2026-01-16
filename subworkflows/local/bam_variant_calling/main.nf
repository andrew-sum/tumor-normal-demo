/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BAM_VARIANT_CALLING SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Somatic variant calling with GATK Mutect2
    - Pairs tumor samples with matched normal by patient ID
    - Runs Mutect2 in tumor-normal mode
    - Filters variants with FilterMutectCalls
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GATK4_MUTECT2            } from '../../../modules/nf-core/gatk4/mutect2/main'
include { GATK4_FILTERMUTECTCALLS  } from '../../../modules/nf-core/gatk4/filtermutectcalls/main'

workflow BAM_VARIANT_CALLING {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta            // channel: [ val(meta), path(fasta) ]
    ch_fai              // channel: [ val(meta), path(fai) ]
    ch_dict             // channel: [ val(meta), path(dict) ]
    ch_germline_resource     // channel: [ path(vcf) ]
    ch_germline_resource_tbi // channel: [ path(tbi) ]
    ch_panel_of_normals      // channel: [ path(vcf) ]
    ch_panel_of_normals_tbi  // channel: [ path(tbi) ]
    ch_intervals        // channel: [ path(bed) ]

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
    // Output: [ meta_tumor, tumor_bam, tumor_bai, normal_bam, normal_bai, normal_sample_name ]
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
            // Create new meta with both sample names for Mutect2
            def new_meta = tumor_meta.clone()
            new_meta.normal_id = normal_meta.id
            [ new_meta, [ tumor_bam, normal_bam ], [ tumor_bai, normal_bai ] ]
        }
        .set { ch_paired }

    //
    // Prepare input for Mutect2
    // Input format: [ meta, [bams], [bais], intervals ]
    //
    ch_paired
        .combine(ch_intervals.map { it -> [ it ] }.ifEmpty( [ [] ] ))
        .map { meta, bams, bais, intervals ->
            [ meta, bams, bais, intervals ]
        }
        .set { ch_mutect2_input }

    //
    // Prepare fai channel - module expects [ meta, fai, gzi ]
    // gzi is only needed for bgzipped fasta, pass empty for regular fasta
    //
    ch_fai_gzi = ch_fai.map { meta, fai -> [ meta, fai, [] ] }

    //
    // MODULE: GATK4_MUTECT2
    //
    GATK4_MUTECT2 (
        ch_mutect2_input,
        ch_fasta,
        ch_fai_gzi,
        ch_dict,
        [],  // alleles
        [],  // alleles_tbi
        ch_germline_resource.ifEmpty( [] ),
        ch_germline_resource_tbi.ifEmpty( [] ),
        ch_panel_of_normals.ifEmpty( [] ),
        ch_panel_of_normals_tbi.ifEmpty( [] )
    )
    ch_versions = ch_versions.mix(GATK4_MUTECT2.out.versions.first())

    //
    // Prepare input for FilterMutectCalls
    // Input format: [ meta, vcf, tbi, stats, orientationbias, segmentation, table, estimate ]
    //
    GATK4_MUTECT2.out.vcf
        .join(GATK4_MUTECT2.out.tbi)
        .join(GATK4_MUTECT2.out.stats)
        .map { meta, vcf, tbi, stats ->
            // Pass empty values for optional contamination inputs
            [ meta, vcf, tbi, stats, [], [], [], [] ]
        }
        .set { ch_filter_input }

    //
    // Prepare fai for filter - just needs [ meta, fai ]
    //
    ch_fai_simple = ch_fai.map { meta, fai -> [ meta, fai ] }

    //
    // MODULE: GATK4_FILTERMUTECTCALLS
    //
    GATK4_FILTERMUTECTCALLS (
        ch_filter_input,
        ch_fasta,
        ch_fai_simple,
        ch_dict
    )
    ch_versions = ch_versions.mix(GATK4_FILTERMUTECTCALLS.out.versions.first())

    emit:
    vcf_unfiltered = GATK4_MUTECT2.out.vcf                  // channel: [ val(meta), path(vcf) ]
    vcf            = GATK4_FILTERMUTECTCALLS.out.vcf        // channel: [ val(meta), path(vcf) ]
    tbi            = GATK4_FILTERMUTECTCALLS.out.tbi        // channel: [ val(meta), path(tbi) ]
    stats          = GATK4_FILTERMUTECTCALLS.out.stats      // channel: [ val(meta), path(stats) ]
    versions       = ch_versions                            // channel: [ path(versions.yml) ]
}
