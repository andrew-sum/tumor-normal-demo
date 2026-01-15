# nf-core Pipeline Development Guide

## Directory Structure

```
modules/nf-core/<tool>/main.nf      # Tool modules
modules/local/<tool>/main.nf        # Custom modules
subworkflows/nf-core/<name>/main.nf # Reusable subworkflows
subworkflows/local/<name>/main.nf   # Custom subworkflows
workflows/<pipeline>.nf             # Main workflow
conf/modules.config                 # Module arguments and publishing
```

## Module Structure (DSL2)

```nextflow
process TOOL_NAME {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/<container>' :
        'biocontainers/<container>' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tool_command \\
        $args \\
        --input $reads \\
        --output ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        toolname: \$(tool --version | sed 's/.*//g')
    END_VERSIONS
    """
}
```

## Process Labels (conf/base.config)

| Label | CPUs | Memory | Time |
|-------|------|--------|------|
| `process_single` | 1 | 6 GB | 4h |
| `process_low` | 2 | 12 GB | 4h |
| `process_medium` | 6 | 36 GB | 8h |
| `process_high` | 12 | 72 GB | 16h |
| `process_high_memory` | 12 | 200 GB | 16h |

## Meta Map Convention

Always pass sample metadata as first element in tuples:
```nextflow
// meta contains at minimum: [ id: 'sample_name', single_end: false ]
tuple val(meta), path(reads)
```

## Channel Operations

```nextflow
// Branch by condition
ch_reads.branch { meta, reads ->
    single: meta.single_end
    paired: !meta.single_end
}

// Join channels by meta.id
ch_bam.join(ch_bai)

// Combine with reference
ch_reads.combine(ch_reference)

// Collect files for MultiQC
ch_reports.collect{ it[1] }.ifEmpty([])
```

## Module Config (conf/modules.config)

```nextflow
process {
    withName: 'FASTQC' {
        ext.args = '--quiet'
        publishDir = [
            path: { "${params.outdir}/fastqc" },
            mode: params.publish_dir_mode,
            saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
        ]
    }
}
```

## Workflow Pattern

```nextflow
include { TOOL_A } from '../modules/nf-core/tool_a/main'
include { TOOL_B } from '../modules/nf-core/tool_b/main'

workflow PIPELINE {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    TOOL_A(ch_samplesheet)
    ch_versions = ch_versions.mix(TOOL_A.out.versions.first())

    TOOL_B(TOOL_A.out.bam)
    ch_versions = ch_versions.mix(TOOL_B.out.versions.first())

    emit:
    bam      = TOOL_B.out.bam
    versions = ch_versions
}
```

## Subworkflow Pattern

Subworkflows chain multiple modules into reusable units:

```nextflow
// subworkflows/local/bam_sort_index/main.nf
include { SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_SORT_INDEX {
    take:
    ch_bam  // channel: [ val(meta), path(bam) ]

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_SORT(ch_bam)
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = SAMTOOLS_SORT.out.bam    // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai   // channel: [ val(meta), path(bai) ]
    versions = ch_versions              // channel: [ path(versions.yml) ]
}
```

## Using Subworkflows

```nextflow
// In main workflow
include { BAM_SORT_INDEX } from '../subworkflows/local/bam_sort_index/main'

workflow PIPELINE {
    take:
    ch_samplesheet

    main:
    ch_versions = Channel.empty()

    ALIGN(ch_samplesheet)
    ch_versions = ch_versions.mix(ALIGN.out.versions)

    BAM_SORT_INDEX(ALIGN.out.bam)
    ch_versions = ch_versions.mix(BAM_SORT_INDEX.out.versions)

    // Join bam and bai for downstream tools
    ch_bam_bai = BAM_SORT_INDEX.out.bam.join(BAM_SORT_INDEX.out.bai)

    emit:
    bam_bai  = ch_bam_bai
    versions = ch_versions
}
```

## Subworkflow with Optional Inputs

```nextflow
workflow ALIGN_AND_QC {
    take:
    ch_reads      // required: [ val(meta), path(reads) ]
    ch_reference  // required: [ path(fasta), path(index) ]
    skip_qc       // val: boolean

    main:
    ch_versions = Channel.empty()
    ch_reports  = Channel.empty()

    ALIGNER(ch_reads, ch_reference)
    ch_versions = ch_versions.mix(ALIGNER.out.versions.first())

    if (!skip_qc) {
        SAMTOOLS_STATS(ALIGNER.out.bam)
        ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.first())
        ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats)
    }

    emit:
    bam      = ALIGNER.out.bam
    reports  = ch_reports
    versions = ch_versions
}
```

## Version Tracking

Collect versions via Nextflow topic:
```nextflow
ch_versions.collectFile(name: 'versions.yml', storeDir: "${params.outdir}/pipeline_info")
```

## Common nf-core Commands

```bash
nf-core modules install <tool>          # Install module from nf-core
nf-core modules list remote             # List available modules
nf-core pipelines lint                   # Check pipeline compliance
nextflow run . -profile test,docker     # Test pipeline
```

## Key Rules

1. Always emit `versions.yml` from every process
2. Use `meta` map for sample tracking through pipeline
3. Use `task.ext.args` for user-configurable arguments
4. Use `task.ext.prefix` for output file naming
5. Quote all variables in shell blocks
6. Use `\$` for shell variables, `$` for Nextflow variables
7. Never hardcode paths - use `params` or `moduleDir`
8. Use micromamba nextflow environment for nextflow and nf-core modules