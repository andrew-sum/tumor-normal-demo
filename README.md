# andrew-sum/tumor-normal-demo

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.5.1-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.5.1)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**andrew-sum/tumor-normal-demo** is a Nextflow bioinformatics pipeline for somatic variant and copy number analysis from matched tumor-normal whole exome or whole genome sequencing data. It takes paired-end FASTQ files for a tumor and matched normal sample, performs read QC, alignment, preprocessing, and calls somatic SNVs/indels and copy number variants.

**Pipeline steps:**

1. Read QC and trimming ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [`fastp`](https://github.com/OpenGEN/fastp))
2. Alignment to reference genome ([`BWA-MEM2`](https://github.com/bwa-mem2/bwa-mem2))
3. Duplicate marking ([`GATK MarkDuplicates`](https://gatk.broadinstitute.org/))
4. Base quality score recalibration — BQSR ([`GATK BaseRecalibrator`](https://gatk.broadinstitute.org/), [`GATK ApplyBQSR`](https://gatk.broadinstitute.org/))
5. Coverage QC ([`mosdepth`](https://github.com/brentp/mosdepth))
6. Somatic SNV/indel calling ([`GATK Mutect2`](https://gatk.broadinstitute.org/), [`GATK FilterMutectCalls`](https://gatk.broadinstitute.org/))
7. Copy number variant calling ([`CNVkit`](https://cnvkit.readthedocs.io/))
8. Aggregated QC report ([`MultiQC`](http://multiqc.info/))

## Usage

Prepare a samplesheet `samplesheet.csv` with your tumor-normal pair:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
PATIENT1,XX,0,PATIENT1_normal,1,normal_R1.fastq.gz,normal_R2.fastq.gz
PATIENT1,XX,1,PATIENT1_tumor,1,tumor_R1.fastq.gz,tumor_R2.fastq.gz
```

- `status`: `0` = normal, `1` = tumor
- Each row is one sample (or lane); multiple lanes for the same sample are merged automatically

Run the pipeline:

```bash
nextflow run andrew-sum/tumor-normal-demo \
   -profile docker \
   --input samplesheet.csv \
   --outdir results \
   --genome GRCh38
```

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Key parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--input` | Samplesheet CSV | required |
| `--outdir` | Output directory | required |
| `--genome` | Reference genome (iGenomes key) | required |
| `--intervals` | Target BED file (exome/panel) | optional |
| `--skip_markduplicates` | Skip duplicate marking + BQSR | false |
| `--skip_variant_calling` | Skip Mutect2 somatic calling | false |
| `--skip_cnv_calling` | Skip CNVkit CNV calling | false |

## Output

```
results/
├── fastqc/             # Raw read QC
├── fastp/              # Trimming reports and trimmed reads
├── alignment/          # Sorted BAM + index
├── preprocessing/
│   ├── markduplicates/ # Duplicate metrics
│   └── recalibrated/   # BQSR BAM + index
├── reports/mosdepth/   # Coverage summaries
├── variant_calling/mutect2/  # Somatic VCFs (raw + filtered)
├── cnv_calling/cnvkit/ # CNV segments, scatter/diagram plots
├── multiqc/            # Aggregated QC report
└── pipeline_info/      # Software versions, execution trace
```

## Credits

andrew-sum/tumor-normal-demo was written by Andrew Sum.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
