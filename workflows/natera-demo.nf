/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQ_TRIM_QC          } from '../subworkflows/local/fastq_trim_qc/main'
include { FASTQ_ALIGN            } from '../subworkflows/local/fastq_align/main'
include { BAM_MARKDUP_BQSR       } from '../subworkflows/local/bam_markdup_bqsr/main'
include { BAM_VARIANT_CALLING    } from '../subworkflows/local/bam_variant_calling/main'
include { BAM_CNV_CALLING        } from '../subworkflows/local/bam_cnv_calling/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_natera-demo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NATERA_DEMO {

    take:
    ch_reads                 // channel: [ val(meta), [ path(reads) ] ]
    ch_fasta                 // channel: [ val(meta), path(fasta) ]
    ch_bwamem2               // channel: [ val(meta), path(index) ]
    ch_fai                   // channel: [ val(meta), path(fai) ]
    ch_dict                  // channel: [ val(meta), path(dict) ]
    ch_known_sites           // channel: [ val(meta), path(vcf) ]
    ch_known_sites_tbi       // channel: [ val(meta), path(tbi) ]
    ch_intervals             // channel: [ path(bed) ]
    ch_germline_resource     // channel: [ path(vcf) ] - optional gnomAD
    ch_germline_resource_tbi // channel: [ path(tbi) ] - optional
    ch_panel_of_normals      // channel: [ path(vcf) ] - optional PoN
    ch_panel_of_normals_tbi  // channel: [ path(tbi) ] - optional

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: QC and trim raw reads
    //
    FASTQ_TRIM_QC ( ch_reads )

    ch_trimmed_reads = FASTQ_TRIM_QC.out.reads
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_QC.out.fastqc_zip.collect{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_TRIM_QC.out.fastp_json.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQ_TRIM_QC.out.versions)

    //
    // SUBWORKFLOW: Align reads to reference genome
    //
    FASTQ_ALIGN (
        ch_trimmed_reads,
        ch_fasta,
        ch_bwamem2
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN.out.versions)

    //
    // SUBWORKFLOW: Mark duplicates and base quality score recalibration
    //
    if (!params.skip_markduplicates) {
        BAM_MARKDUP_BQSR (
            FASTQ_ALIGN.out.bam_bai,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_known_sites,
            ch_known_sites_tbi,
            ch_intervals
        )
        ch_versions = ch_versions.mix(BAM_MARKDUP_BQSR.out.versions)

        // Add markdup metrics and mosdepth to MultiQC
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUP_BQSR.out.metrics.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(BAM_MARKDUP_BQSR.out.coverage_global.collect{it[1]})

        ch_bam_bai = BAM_MARKDUP_BQSR.out.bam_bai
    } else {
        ch_bam_bai = FASTQ_ALIGN.out.bam_bai
    }

    //
    // SUBWORKFLOW: Somatic variant calling with Mutect2
    //
    ch_vcf = Channel.empty()
    ch_vcf_tbi = Channel.empty()

    if (!params.skip_variant_calling) {
        BAM_VARIANT_CALLING (
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_germline_resource,
            ch_germline_resource_tbi,
            ch_panel_of_normals,
            ch_panel_of_normals_tbi,
            ch_intervals
        )
        ch_versions = ch_versions.mix(BAM_VARIANT_CALLING.out.versions)

        ch_vcf = BAM_VARIANT_CALLING.out.vcf
        ch_vcf_tbi = BAM_VARIANT_CALLING.out.tbi
    }

    //
    // SUBWORKFLOW: CNV calling with CNVkit
    //
    ch_cnv_cns = Channel.empty()

    if (!params.skip_cnv_calling) {
        BAM_CNV_CALLING (
            ch_bam_bai,
            ch_fasta,
            ch_fai,
            ch_intervals
        )
        ch_versions = ch_versions.mix(BAM_CNV_CALLING.out.versions)

        ch_cnv_cns = BAM_CNV_CALLING.out.cns
    }

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'natera-demo_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    trimmed_reads  = ch_trimmed_reads            // channel: [ val(meta), [ path(reads) ] ]
    bam_bai        = ch_bam_bai                  // channel: [ val(meta), path(bam), path(bai) ]
    vcf            = ch_vcf                      // channel: [ val(meta), path(vcf) ]
    vcf_tbi        = ch_vcf_tbi                  // channel: [ val(meta), path(tbi) ]
    cnv_cns        = ch_cnv_cns                  // channel: [ val(meta), path(cns) ]
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
