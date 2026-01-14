# Implementation Todos

## 0. Input Schema Setup ✅
- [x] Update `assets/schema_input.json` with sarek-compatible columns:
  - patient, sex, status, sample, lane, fastq_1, fastq_2
- [x] Update samplesheet validation in `subworkflows/local/utils_nfcore_natera-demo_pipeline/`
- [x] Create test samplesheet `assets/samplesheet_test.csv` (HCC1395 data)
- [x] Fix workflow naming (NATERA-DEMO → NATERA_DEMO)
- [x] Test pipeline with FastQC - **PASSED**

---

## 1. FASTQ_TRIM_QC ✅
**Path**: `subworkflows/local/fastq_trim_qc/main.nf`

**Modules to install**:
- [x] (fastqc already installed)
- [x] `nf-core modules install fastp`

**Implementation**:
- [x] Create subworkflow: FASTQC → FASTP
- [x] Configure in `conf/modules.config`
- [x] Test subworkflow independently

---

## 2. FASTQ_ALIGN
**Path**: `subworkflows/local/fastq_align/main.nf`

**Modules to install**:
- [ ] `nf-core modules install bwa/mem2`
- [ ] `nf-core modules install samtools/sort`
- [ ] `nf-core modules install samtools/index`

**Implementation**:
- [ ] Create subworkflow: BWA-MEM2 → SAMTOOLS_SORT → SAMTOOLS_INDEX
- [ ] Configure in `conf/modules.config`
- [ ] Test subworkflow independently

---

## 3. BAM_MARKDUP_BQSR
**Path**: `subworkflows/local/bam_markdup_bqsr/main.nf`

**Modules to install**:
- [ ] `nf-core modules install gatk4/markduplicates`
- [ ] `nf-core modules install gatk4/baserecalibrator`
- [ ] `nf-core modules install gatk4/applybqsr`
- [ ] `nf-core modules install mosdepth`

**Implementation**:
- [ ] Create subworkflow: MARKDUPLICATES → BASERECALIBRATOR → APPLYBQSR → MOSDEPTH
- [ ] Configure in `conf/modules.config`
- [ ] Test subworkflow independently

---

## 4. BAM_VARIANT_CALLING
**Path**: `subworkflows/local/bam_variant_calling/main.nf`

**Modules to install**:
- [ ] `nf-core modules install gatk4/mutect2`
- [ ] `nf-core modules install gatk4/filtermutectcalls`

**Implementation**:
- [ ] Create subworkflow: MUTECT2 → FILTERMUTECTCALLS
- [ ] Handle tumor-normal pairing input
- [ ] Configure in `conf/modules.config`
- [ ] Test subworkflow independently

---

## 5. BAM_CNV_CALLING
**Path**: `subworkflows/local/bam_cnv_calling/main.nf`

**Modules to install**:
- [ ] `nf-core modules install cnvkit/batch`

**Implementation**:
- [ ] Create subworkflow: CNVKIT_BATCH
- [ ] Configure in `conf/modules.config`
- [ ] Test subworkflow independently

---

## 6. BAM_MSI
**Path**: `subworkflows/local/bam_msi/main.nf`

**Modules to install**:
- [ ] `nf-core modules install msisensorpro/scan`
- [ ] `nf-core modules install msisensorpro/msi_somatic`

**Implementation**:
- [ ] Create subworkflow: MSISENSORPRO_SCAN → MSISENSORPRO_MSI_SOMATIC
- [ ] Configure in `conf/modules.config`
- [ ] Test subworkflow independently

---

## 7. VCF_ANNOTATE_VEP
**Path**: `subworkflows/local/vcf_annotate_vep/main.nf`

**Modules to install**:
- [ ] `nf-core modules install ensemblvep/vep`
- [ ] `nf-core modules install bcftools/stats`

**Implementation**:
- [ ] Create subworkflow: VEP → BCFTOOLS_STATS
- [ ] Configure in `conf/modules.config`
- [ ] Test subworkflow independently

---

## 8. PIPELINE_REPORTING
**Path**: `subworkflows/local/pipeline_reporting/main.nf`

**Modules to install**:
- [ ] (multiqc already installed)

**Implementation**:
- [ ] Create subworkflow: Collect all QC files → MULTIQC
- [ ] Aggregate: FastQC, fastp, mosdepth, bcftools stats, Mutect2 stats
- [ ] Configure custom MultiQC config in `assets/multiqc_config.yml`
- [ ] Test subworkflow independently

---

## 9. Main Workflow Integration
**Path**: `workflows/natera-demo.nf`

- [ ] Wire all subworkflows together
- [ ] Add tumor-normal channel pairing logic
- [ ] Update parameters in `nextflow.config`
- [ ] Update `nextflow_schema.json`

---

## 10. Final Testing
- [ ] Run full pipeline with test data
- [ ] Verify all outputs
- [ ] Run `nf-core pipelines lint`
- [ ] Fix linting issues
