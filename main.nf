#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    andrew/natera-demo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/andrew/natera-demo
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { NATERA_DEMO  } from './workflows/natera-demo'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_natera-demo_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_natera-demo_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_natera-demo_pipeline'
include { BWAMEM2_INDEX                  } from './modules/nf-core/bwamem2/index/main'
include { SAMTOOLS_FAIDX                 } from './modules/nf-core/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY } from './modules/nf-core/gatk4/createsequencedictionary/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Reference genome files from igenomes
params.fasta   = getGenomeAttribute('fasta')
params.bwamem2 = getGenomeAttribute('bwamem2')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow ANDREW_NATERA_DEMO {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:

    //
    // Set up reference channels
    //
    ch_fasta = params.fasta ? Channel.fromPath(params.fasta).map { [[id: 'reference'], it] }.collect() : Channel.empty()

    //
    // Build BWA-MEM2 index if not provided
    //
    if (params.bwamem2) {
        ch_bwamem2 = Channel.fromPath(params.bwamem2).map { [[id: 'reference'], it] }.collect()
    } else {
        BWAMEM2_INDEX(ch_fasta)
        ch_bwamem2 = BWAMEM2_INDEX.out.index
    }

    //
    // Generate FASTA index (.fai) if not provided
    //
    if (params.fasta_fai) {
        ch_fai = Channel.fromPath(params.fasta_fai).map { [[id: 'reference'], it] }.collect()
    } else {
        SAMTOOLS_FAIDX(ch_fasta, [[],[]], false)
        ch_fai = SAMTOOLS_FAIDX.out.fai
    }

    //
    // Generate sequence dictionary (.dict) if not provided
    //
    if (params.dict) {
        ch_dict = Channel.fromPath(params.dict).map { [[id: 'reference'], it] }.collect()
    } else {
        GATK4_CREATESEQUENCEDICTIONARY(ch_fasta)
        ch_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict
    }

    //
    // Set up known sites for BQSR (optional)
    //
    ch_known_sites = params.known_sites
        ? Channel.fromPath(params.known_sites).map { [[id: 'known_sites'], it] }.collect()
        : Channel.value([[],[]])
    ch_known_sites_tbi = params.known_sites_tbi
        ? Channel.fromPath(params.known_sites_tbi).map { [[id: 'known_sites'], it] }.collect()
        : Channel.value([[],[]])

    //
    // Set up intervals (optional)
    //
    ch_intervals = params.intervals
        ? Channel.fromPath(params.intervals).collect()
        : Channel.value([])

    //
    // WORKFLOW: Run pipeline
    //
    NATERA_DEMO (
        samplesheet,
        ch_fasta,
        ch_bwamem2,
        ch_fai,
        ch_dict,
        ch_known_sites,
        ch_known_sites_tbi,
        ch_intervals
    )
    emit:
    multiqc_report = NATERA_DEMO.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    ANDREW_NATERA_DEMO (
        PIPELINE_INITIALISATION.out.samplesheet
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
        ANDREW_NATERA_DEMO.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
