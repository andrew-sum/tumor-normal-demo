//
// BAM_MARKDUP_BQSR: Mark duplicates and perform base quality score recalibration
//

include { GATK4_MARKDUPLICATES   } from '../../../modules/nf-core/gatk4/markduplicates/main'
include { GATK4_BASERECALIBRATOR } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_APPLYBQSR        } from '../../../modules/nf-core/gatk4/applybqsr/main'
include { MOSDEPTH               } from '../../../modules/nf-core/mosdepth/main'

workflow BAM_MARKDUP_BQSR {

    take:
    ch_bam_bai         // channel: [ val(meta), path(bam), path(bai) ]
    ch_fasta           // channel: [ val(meta), path(fasta) ]
    ch_fai             // channel: [ val(meta), path(fai) ]
    ch_dict            // channel: [ val(meta), path(dict) ]
    ch_known_sites     // channel: [ val(meta), path(vcf) ]
    ch_known_sites_tbi // channel: [ val(meta), path(tbi) ]
    ch_intervals       // channel: [ path(bed) ] - optional target regions

    main:
    ch_versions = Channel.empty()
    ch_table    = Channel.empty()

    //
    // Prepare inputs - extract just the BAM for markduplicates
    //
    ch_bam = ch_bam_bai.map { meta, bam, bai -> [ meta, bam ] }

    //
    // MODULE: Mark duplicates
    //
    GATK4_MARKDUPLICATES (
        ch_bam,
        ch_fasta.map { meta, fasta -> fasta },
        ch_fai.map { meta, fai -> fai }
    )
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions.first())

    //
    // Prepare marked BAM + BAI
    //
    ch_marked_bam_bai = GATK4_MARKDUPLICATES.out.bam
        .join(GATK4_MARKDUPLICATES.out.bai)

    ch_marked_bam_bai.view()

    //
    // BQSR: Only run if not skipped via params
    //
    if (!params.skip_bqsr) {
        //
        // Prepare input for BaseRecalibrator: [ meta, bam, bai, intervals ]
        //
        ch_bam_for_bqsr = ch_marked_bam_bai
            .map { meta, bam, bai ->
                [ meta, bam, bai, [] ]  // empty intervals
            }

        //
        // MODULE: Build recalibration table
        //
        GATK4_BASERECALIBRATOR (
            ch_bam_for_bqsr,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_known_sites,
            ch_known_sites_tbi
        )
        ch_versions = ch_versions.mix(GATK4_BASERECALIBRATOR.out.versions.first())
        ch_table = GATK4_BASERECALIBRATOR.out.table

        //
        // Prepare input for ApplyBQSR: [ meta, bam, bai, table, intervals ]
        //
        ch_bam_for_apply = ch_marked_bam_bai
            .join(GATK4_BASERECALIBRATOR.out.table)
            .map { meta, bam, bai, table ->
                [ meta, bam, bai, table, [] ]
            }

        //
        // MODULE: Apply recalibration
        //
        GATK4_APPLYBQSR (
            ch_bam_for_apply,
            ch_fasta.map { meta, fasta -> fasta },
            ch_fai.map { meta, fai -> fai },
            ch_dict.map { meta, dict -> dict }
        )
        ch_versions = ch_versions.mix(GATK4_APPLYBQSR.out.versions.first())

        // Use recalibrated BAMs
        ch_final_bam_bai = GATK4_APPLYBQSR.out.bam
            .join(GATK4_APPLYBQSR.out.bai)
    } else {
        // Skip BQSR - use marked duplicates BAMs
        ch_final_bam_bai = ch_marked_bam_bai
    }

    //
    // Prepare input for Mosdepth: [ meta, bam, bai, bed ]
    //
    ch_bam_for_mosdepth = ch_final_bam_bai
        .map { meta, bam, bai ->
            [ meta, bam, bai, [] ]
        }

    //
    // MODULE: Calculate coverage statistics
    //
    MOSDEPTH (
        ch_bam_for_mosdepth,
        ch_fasta
    )

    emit:
    bam              = ch_final_bam_bai.map { meta, bam, bai -> [meta, bam] }
    bai              = ch_final_bam_bai.map { meta, bam, bai -> [meta, bai] }
    bam_bai          = ch_final_bam_bai
    metrics          = GATK4_MARKDUPLICATES.out.metrics
    table            = ch_table
    coverage_summary = MOSDEPTH.out.summary_txt
    coverage_global  = MOSDEPTH.out.global_txt
    versions         = ch_versions
}
