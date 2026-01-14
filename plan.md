# Exome Variant Calling Pipeline Plan

## Overview

Build a tumor-normal paired exome sequencing pipeline for detecting somatic mutations (SNVs, Indels, CNVs) with clinical-grade annotation, similar to Natera's Altera workflow.

## References

### nf-core Reference Pipeline

**nf-core/sarek** - The gold-standard nf-core pipeline for somatic/germline variant calling
- GitHub: https://github.com/nf-core/sarek
- Documentation: https://nf-co.re/sarek
- Our pipeline follows sarek's architecture for tumor-normal paired analysis

### Literature

| Reference | Description |
|-----------|-------------|
| [Garcia M, et al. F1000Research 2020](https://f1000research.com/articles/9-63) | Sarek pipeline paper - "Sarek: A portable workflow for whole-genome sequencing analysis of germline and somatic variants" |
| [Cibulskis K, et al. Nature Biotechnology 2013](https://www.nature.com/articles/nbt.2514) | Original MuTect paper - foundational somatic variant calling methodology |
| [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731-Somatic-short-variant-discovery-SNVs-Indels) | Somatic short variant discovery workflow documentation |
| [Fang LT, et al. Nature Biotechnology 2021](https://www.nature.com/articles/s41587-021-00993-6) | SEQC2 somatic mutation benchmarking - established truth sets for validation |

## Test Datasets

### SEQC2 Benchmark Data (HCC1395)

The FDA-led SEQC2 consortium provides the gold-standard tumor-normal benchmark. We use the same data as nf-core/sarek.

| Sample | Type | Cell Line | SRA |
|--------|------|-----------|-----|
| **HCC1395T** | Tumor | HCC1395 (Triple-negative breast cancer) | SRR7890918 |
| **HCC1395N** | Normal | HCC1395BL (B-lymphocyte) | SRR7890919 |

**Ground truth**: 39,536 SNVs + 2,020 indels with high-confidence regions

**Test samplesheet** (`assets/samplesheet_test.csv`):
```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
HCC1395,XX,0,HCC1395N,1,s3://ngi-igenomes/test-data/sarek/SRR7890919_WES_HCC1395BL-EA_normal_1.fastq.gz,s3://ngi-igenomes/test-data/sarek/SRR7890919_WES_HCC1395BL-EA_normal_2.fastq.gz
HCC1395,XX,1,HCC1395T,1,s3://ngi-igenomes/test-data/sarek/SRR7890918_WES_HCC1395-EA_tumor_1.fastq.gz,s3://ngi-igenomes/test-data/sarek/SRR7890918_WES_HCC1395-EA_tumor_2.fastq.gz
```

**Data characteristics**:
- WES: ~100x tumor, ~80x normal coverage
- Hosted on AWS S3 (ngi-igenomes bucket)
- Same data used by nf-core/sarek for validation

### Quick CI Testing

For rapid development testing, use nf-core's downsampled datasets:
```bash
nextflow run . -profile test,docker --outdir results
```

### Alternative: Direct SRA Download

If S3 access is unavailable:
```bash
# Using sra-tools
prefetch SRR7890918 SRR7890919
fasterq-dump SRR7890918 SRR7890919 --split-files
```

## Comparison with nf-core/sarek

| Aspect | Our Pipeline | nf-core/sarek |
|--------|--------------|---------------|
| **Scope** | Tumor-normal WES focused | WGS, WES, targeted; germline + somatic |
| **Complexity** | ~15 modules, easy to customize | 50+ modules, full-featured |
| **Use case** | Demo/learning/focused analysis | Production clinical pipelines |

### Feature Comparison

| Feature | Ours | Sarek |
|---------|------|-------|
| **Aligners** | BWA-MEM2 | BWA-MEM, BWA-MEM2, dragmap, Sentieon |
| **SNV/Indel** | Mutect2 | Mutect2, Strelka2, Freebayes, DeepVariant, MuSE |
| **CNV** | CNVkit | CNVkit, ASCAT, Control-FREEC |
| **SV** | ❌ | Manta, TIDDIT |
| **MSI** | ✅ MSIsensor-pro | ✅ MSIsensor-pro |
| **Coverage** | ✅ mosdepth | ✅ mosdepth |
| **UMI** | ❌ | ✅ fgbio |
| **Annotation** | VEP | VEP, SnpEff |
| **Germline** | ❌ | ✅ HaplotypeCaller, DeepVariant |

### Why Build Our Own?

1. **Simplicity** - Easier to understand, modify, and debug
2. **Demo-friendly** - Clear workflow for interviews/presentations
3. **Learning** - Understand nf-core patterns without complexity
4. **Customization** - Add custom steps without navigating large codebase

## Pipeline Architecture

```
FASTQ (Tumor + Normal)
    │
    ▼
┌─────────────────────────────────────┐
│  1. FASTQ_TRIM_QC                   │
│     ├── FASTQC                      │
│     └── FASTP                       │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  2. FASTQ_ALIGN                     │
│     ├── BWA-MEM2                    │
│     ├── SAMTOOLS SORT               │
│     └── SAMTOOLS INDEX              │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  3. BAM_MARKDUP_BQSR                │
│     ├── GATK MARKDUPLICATES         │
│     ├── GATK BASERECALIBRATOR       │
│     ├── GATK APPLYBQSR              │
│     └── MOSDEPTH                    │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  4. BAM_VARIANT_CALLING             │
│     ├── GATK MUTECT2                │
│     └── GATK FILTERMUTECTCALLS      │
├─────────────────────────────────────┤
│  5. BAM_CNV_CALLING                 │
│     └── CNVKIT                      │
├─────────────────────────────────────┤
│  6. BAM_MSI                         │
│     ├── MSISENSORPRO SCAN           │
│     └── MSISENSORPRO MSI_SOMATIC    │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  7. VCF_ANNOTATE_VEP                │
│     ├── ENSEMBL VEP                 │
│     └── BCFTOOLS STATS              │
└─────────────────────────────────────┘
    │
    ▼
┌─────────────────────────────────────┐
│  8. PIPELINE_REPORTING              │
│     └── MULTIQC                     │
└─────────────────────────────────────┘
```

## Input Samplesheet Format

Using nf-core/sarek compatible format:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
HCC1395,XX,0,HCC1395N,1,normal_R1.fastq.gz,normal_R2.fastq.gz
HCC1395,XX,1,HCC1395T,1,tumor_R1.fastq.gz,tumor_R2.fastq.gz
```

| Column | Description |
|--------|-------------|
| `patient` | Patient/case ID - groups tumor-normal pairs |
| `sex` | XX or XY (for CNV calling) |
| `status` | 0 = normal, 1 = tumor |
| `sample` | Unique sample name |
| `lane` | Sequencing lane (for multi-lane samples) |
| `fastq_1` | Path to R1 FASTQ |
| `fastq_2` | Path to R2 FASTQ |

## Implementation Steps

### Step 1: Update Input Schema
- **File**: `assets/schema_input.json`
- Add columns for tumor/normal paired FASTQs
- Update samplesheet validation

### Step 2: Create Subworkflows

#### 2a: FASTQ_TRIM_QC (`subworkflows/local/fastq_trim_qc/`)
- FASTQC → FASTP
- Takes: reads channel
- Emits: trimmed reads, QC reports, versions

#### 2b: FASTQ_ALIGN (`subworkflows/local/fastq_align/`)
- BWA-MEM2 → SAMTOOLS_SORT → SAMTOOLS_INDEX
- Takes: trimmed reads, reference index
- Emits: sorted BAM + BAI, versions

#### 2c: BAM_MARKDUP_BQSR (`subworkflows/local/bam_markdup_bqsr/`)
- GATK_MARKDUPLICATES → GATK_BASERECALIBRATOR → GATK_APPLYBQSR → MOSDEPTH
- Takes: sorted BAM, reference, known sites
- Emits: analysis-ready BAM, coverage stats, versions

#### 2d: BAM_VARIANT_CALLING (`subworkflows/local/bam_variant_calling/`)
- GATK_MUTECT2 → GATK_FILTERMUTECTCALLS
- Takes: tumor BAM, normal BAM, reference, intervals
- Emits: filtered VCF, versions

#### 2e: BAM_CNV_CALLING (`subworkflows/local/bam_cnv_calling/`)
- CNVKIT_BATCH
- Takes: tumor BAM, normal BAM, reference, targets
- Emits: CNV calls, versions

#### 2f: BAM_MSI (`subworkflows/local/bam_msi/`)
- MSISENSORPRO_SCAN → MSISENSORPRO_MSI_SOMATIC
- Takes: tumor BAM, normal BAM, reference
- Emits: MSI score, MSI classification (MSS/MSI-L/MSI-H), versions

#### 2g: VCF_ANNOTATE_VEP (`subworkflows/local/vcf_annotate_vep/`)
- ENSEMBL_VEP → BCFTOOLS_STATS
- Takes: VCF, cache, plugins (ClinVar, COSMIC, gnomAD)
- Emits: annotated VCF, VCF stats, versions

#### 2h: PIPELINE_REPORTING (`subworkflows/local/pipeline_reporting/`)
- MULTIQC (aggregate all QC)
- Takes: FastQC, fastp, mosdepth, bcftools stats outputs
- Emits: MultiQC HTML report, versions

### Step 3: Install nf-core Modules

```bash
# Preprocessing
nf-core modules install fastp
nf-core modules install bwa/mem2
nf-core modules install samtools/sort
nf-core modules install samtools/index

# BAM processing
nf-core modules install gatk4/markduplicates
nf-core modules install gatk4/baserecalibrator
nf-core modules install gatk4/applybqsr
nf-core modules install mosdepth

# Variant calling
nf-core modules install gatk4/mutect2
nf-core modules install gatk4/filtermutectcalls
nf-core modules install cnvkit/batch
nf-core modules install msisensorpro/scan
nf-core modules install msisensorpro/msi_somatic

# Annotation & reporting
nf-core modules install ensemblvep/vep
nf-core modules install bcftools/stats
```

### Step 4: Update Main Workflow
- **File**: `workflows/natera-demo.nf`
- Wire subworkflows together
- Handle tumor-normal pairing logic
- Collect all QC outputs for MultiQC

### Step 5: Add Parameters
- **File**: `nextflow.config` and `nextflow_schema.json`
- `--fasta`: Reference genome
- `--bwa_index`: BWA index path
- `--known_sites`: Known variant sites for BQSR
- `--intervals`: Target regions BED file
- `--vep_cache`: VEP cache directory
- `--skip_cnv`: Optional skip CNV calling
- `--skip_msi`: Optional skip MSI detection

### Step 6: Configure Modules
- **File**: `conf/modules.config`
- Set publishDir for each module
- Configure ext.args for tools

## Reference Files Required

| Parameter | Description |
|-----------|-------------|
| `fasta` | Reference genome (hg38) |
| `bwa_index` | BWA-MEM2 index |
| `known_sites` | dbSNP + Mills indels for BQSR |
| `intervals` | Exome capture targets BED |
| `vep_cache` | VEP cache with ClinVar, COSMIC |

## Verification

1. Run with test profile:
   ```bash
   nextflow run . -profile test,docker --outdir results
   ```

2. Check outputs:
   - `results/fastp/` - Trimming reports
   - `results/bwa/` - Aligned BAMs
   - `results/mosdepth/` - Coverage statistics
   - `results/mutect2/` - Somatic VCFs
   - `results/cnvkit/` - CNV calls
   - `results/msisensorpro/` - MSI scores and classification
   - `results/vep/` - Annotated VCFs
   - `results/multiqc/` - Aggregated report

3. Validate VCF format:
   ```bash
   bcftools view results/vep/*.vcf.gz | head
   ```
