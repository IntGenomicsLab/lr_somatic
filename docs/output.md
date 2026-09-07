# IntGenomicsLab/lrsomatic: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

### Output Example

The pipeline produces per-sample output directories. Two modes exist depending on whether a matched normal sample is provided.

**Tumor-only sample** (no matched normal, `-TO` variant callers):

```
├── Sample ID
│    ├── ascat
│    ├── bamfiles
│    ├── methylation
│    │    └── tumor
│    │        └── modkit_pileup
│    ├── qc
│    │    ├── tumor
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam_rep1
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── nanoplot_aln
│    │    │   ├── nanoplot_ubam_rep1
│    │    │   └── samtools
│    │    └── whatshap_stats
│    ├── variants
│    │   ├── clairsto
│    │   ├── deepsomatic
│    │   ├── deepvariant
│    │   ├── phased
│    │   └── severus
│    ├── vep
│    │   ├── somatic
│    │   └── SVs
│    └── wakhan
```

**Paired tumor + normal sample**:

```
├── Sample ID
│    ├── ascat
│    ├── bamfiles
│    ├── methylation
│    │    ├── tumor
│    │    │   └── modkit_pileup
│    │    └── normal
│    │        └── modkit_pileup
│    ├── qc
│    │    ├── tumor
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam_rep1
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── nanoplot_aln
│    │    │   ├── nanoplot_ubam_rep1
│    │    │   └── samtools
│    │    ├── normal
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam_rep1
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── nanoplot_aln
│    │    │   ├── nanoplot_ubam_rep1
│    │    │   └── samtools
│    │    └── whatshap_stats
│    ├── variants
│    │   ├── clair3
│    │   ├── clairs
│    │   ├── deepsomatic
│    │   ├── deepvariant
│    │   ├── phased
│    │   └── severus
│    ├── vep
│    │   ├── germline
│    │   ├── somatic
│    │   └── SVs
│    └── wakhan
├── pipeline_info
└── multiqc
```

### `ascat`

<details markdown="1">
<summary>Output files</summary>

```
├── ascat
│   ├── sample.before_correction.sample.tumour.germline.png
│   ├── sample.before_correction.sample.tumour.tumour.png
│   ├── sample.after_correction.sample.tumour.germline.png
│   ├── sample.after_correction.sample.tumour.tumour.png
│   ├── sample.cnvs.txt
│   ├── sample.metrics.txt
│   ├── sample.normal_alleleFrequencies_chr(1-22,X).txt
│   ├── sample.purityploidy.txt
│   ├── sample.segments.txt
│   ├── sample.segments_raw.txt
│   ├── sample.tumour_alleleFrequencies_chr(1-22,X).txt
│   ├── sample.tumour_normalBAF_rawBAF.txt
│   ├── sample.tumour_normalBAF.txt
│   ├── sample.tumor_tumourLogR.txt
│   ├── sample.tumour.ASCATprofile.png
│   ├── sample.tumour.ASPCF.png
│   ├── sample.tumour.rawprofile.png
│   ├── sample.tumour.sunrise.png
```

| File                                                  | Description                                                                                                      |
| ----------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `sample.before_correction.sample.tumour.germline.png` | LogR and BAF plots from the normal sample before correction                                                      |
| `sample.before_correction.sample.tumour.tumour.png`   | LogR and BAF plots from the tumor sample before correction                                                       |
| `sample.cnvs.txt`                                     | a tsv file describing each chromosome segment with a copy number alteration and it's major and minor copy number |
| `sample.metrics.txt`                                  | a tsv file describing summary statistics for the sample                                                          |
| `sample.normal_alleleFrequencies_chr(1-22,X).txt`     | a tsv file describing the snp counts for the normal sample at each position and their respective depths          |
| `sample.purityploidy.txt`                             | a tsv file describing the purity and ploidy values of the sample                                                 |
| `sample.segments.txt`                                 | a tsv file describing each chromosome segment and it's major and minor copy number                               |
| `sample.segments_raw.txt`                             | a tsv file describing each chromosome segment and it's major and minor rounded and raw copy number               |
| `sample.tumour_alleleFrequencies_chr(1-22,X).txt`     | a tsv file describing the snp counts for the tumor sample at each position and their respective depths           |
| `sample.tumour_normalBAF_rawBAF.txt`                  | a tsv file with the raw BAF values in the normal sample                                                          |
| `sample.tumour_normalBAF.txt`                         | a tsv file with the BAF values in the normal sample                                                              |
| `sample.tumour_normalLogR.txt`                        | a tsv file with the LogR values in the normal sample                                                             |
| `sample.tumour_tumourBAF_rawBAF.txt`                  | a tsv file with the raw BAF values in the tumor sample                                                           |
| `sample.tumour_tumourBAF.txt`                         | a tsv file with the corrected BAF values in the tumor sample                                                     |
| `sample.tumour_tumourLogR.txt`                        | a tsv file with the corrected LogR values in the tumor sample                                                    |
| `sample.tumour.ASCATprofile.png`                      | a png file with the corrected overall copy number profile with ploidy, purity, and goodness of fit metrics       |
| `sample.tumour.ASPCF.png`                             | a png file with the corrected LogR and BAF plots of the tumor sample                                             |
| `sample.tumour.rawprofile.png`                        | a png file with the raw overall copy number profile with ploidy, purity, and goodness of fit metrics             |
| `sample.tumour.sunrise.png`                           | a png file with a purity and ploidy fit                                                                          |

</details>

### `bamfiles`

<details markdown="1">
<summary>Output files</summary>

```
├── bamfiles
│   ├── sample_normal.bam
│   ├── sample_normal.bam.bai
│   ├── sample_tumor.bam
│   ├── sample_tumor.bam.bai
```

| File                    | Description                                                                                          |
| ----------------------- | ---------------------------------------------------------------------------------------------------- |
| `sample_normal.bam`     | Aligned and haplotagged bam file (with methylation and nucleosome predictions) for the normal sample |
| `sample_normal.bam.bai` | index file for the normal bam file                                                                   |
| ` sample_tumor.bam`     | Aligned and haplotagged bam file (with methylation and nucleosome predictions) for the tumor sample  |
| `sample_tumor.bam.bai`  | index file for the tumor bam file                                                                    |

</details>

### `qc`

<details markdown="1">
<summary>Output files</summary>

QC outputs are placed under `tumor/` for all samples, and additionally under `normal/` for paired tumor + normal samples. `whatshap_stats/` appears at the top level of `qc/`.

```
├── qc
│   ├── tumor
│   │   ├── cramino_aln
│   │   │   ├── sample_tumor_cramino.txt
│   │   ├── cramino_ubam_rep1
│   │   │   ├── sample_tumor_cramino.txt
│   │   ├── fibertoolsrs
│   │   │   ├── sample_qc.txt
│   │   ├── mosdepth
│   │   │   ├── sample.mosdepth.global.dist.txt
│   │   │   ├── sample.mosdepth.summary.txt
│   │   ├── nanoplot_aln
│   │   │   ├── sample_tumor_aln_NanoStats.txt
│   │   │   ├── sample_tumor_aln_NanoPlot-report.html
│   │   ├── nanoplot_ubam_rep1
│   │   │   ├── sample_tumor_ubam_NanoStats.txt
│   │   │   ├── sample_tumor_ubam_NanoPlot-report.html
│   │   ├── samtools
│   │   │   ├── sample.flagstat
│   │   │   ├── sample.idxstats
│   │   │   ├── sample.stats
│   ├── normal                          # paired samples only
│   │   └── [same subdirectories as tumor]
│   ├── whatshap_stats
│   │   ├── sample.stats.tsv
│   │   ├── sample.blocklist.tsv
```

| File                                                         | Description                                                                                                              |
| ------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------ |
| `cramino_aln/sample_{type}_cramino.txt`                      | cramino QC summary statistics for the aligned bam file                                                                   |
| `cramino_ubam_rep1/sample_{type}_cramino.txt`                | cramino QC summary statistics for the unaligned bam files                                                                |
| `fibertoolsrs/sample_qc.txt`                                 | fibertools QC summary for the bam file                                                                                   |
| `mosdepth/sample.mosdepth.global.dist.txt`                   | a cumulative distribution indicating the proportion of total bases that were covered for at least a given coverage value |
| `mosdepth/sample.mosdepth.summary.txt`                       | overall summary file from mosdepth tool                                                                                  |
| `nanoplot_aln/sample_{type}_aln_NanoStats.txt`               | NanoPlot summary statistics for the aligned BAM file                                                                     |
| `nanoplot_aln/sample_{type}_aln_NanoPlot-report.html`        | NanoPlot interactive HTML report for the aligned BAM file                                                                |
| `nanoplot_ubam_rep1/sample_{type}_ubam_NanoStats.txt`        | NanoPlot summary statistics for the unaligned BAM file                                                                   |
| `nanoplot_ubam_rep1/sample_{type}_ubam_NanoPlot-report.html` | NanoPlot interactive HTML report for the unaligned BAM file                                                              |
| `samtools/sample.flagstat`                                   | a summary of the counts of different samtools flags                                                                      |
| `samtools/sample.idxstats`                                   | a summary of the number of mapped and unmapped reads                                                                     |
| `samtools/sample.stats`                                      | summary statistics from the bamfile                                                                                      |
| `whatshap_stats/sample.stats.tsv`                            | WhatsHap phasing statistics per chromosome including phase block N50 and switch error rates                              |
| `whatshap_stats/sample.blocklist.tsv`                        | list of all phase blocks with their genomic coordinates                                                                  |

</details>

### `methylation`

<details markdown="1">
<summary>Output files</summary>

```
├── methylation
│   ├── tumor
│   │   └── modkit_pileup
│   │       └── sample.bed.gz
│   ├── normal                          # paired samples only
│   │   └── modkit_pileup
│   │       └── sample.bed.gz
```

| File                                         | Description                                                                         |
| -------------------------------------------- | ----------------------------------------------------------------------------------- |
| `{tumor,normal}/modkit_pileup/sample.bed.gz` | Modkit pileup BED file containing per-CpG methylation frequency and coverage values |

</details>

### `variants`

<details markdown="1">
<summary>Output files</summary>

#### `clair3`

```
├── clair3
│   ├── merge_output.vcf.gz
│   ├── merge_output.vcf.gz.tbi
```

| File                  | Description                                       |
| --------------------- | ------------------------------------------------- |
| `merge_output.vcf.gz` | Merged germline indel and snv calls in vcf format |
| `merge_output.vcf.gz` | index for germline small variant calls            |

#### `clairS`

Present in **paired** (tumor + normal) samples.

```
├── clairs
│   ├── indel.vcf.gz
│   ├── indel.vcf.gz.tbi
│   ├── snv.vcf.gz
│   ├── snv.vcf.gz.tbi
```

| File               | Description                       |
| ------------------ | --------------------------------- |
| `indel.vcf.gz`     | Somatic indel calls in vcf format |
| `indel.vcf.gz.tbi` | Index for somatic indel calls     |
| `snv.vcf.gz`       | Somatic SNV calls in vcf format   |
| `snv.vcf.gz.tbi`   | Index for somatic SNV calls       |

#### `clairS-TO`

Present in **tumor-only** samples (no matched normal).

```
├── clairsto
│   ├── germline.vcf.gz
│   ├── germline.vcf.gz.tbi
│   ├── indel.vcf.gz
│   ├── indel.vcf.gz.tbi
│   ├── snv.vcf.gz
│   ├── snv.vcf.gz.tbi
│   ├── somatic.vcf.gz
│   ├── somatic.vcf.gz.tbi
```

| File                  | Description                                                           |
| --------------------- | --------------------------------------------------------------------- |
| `germline.vcf.gz`     | SNV and indel calls marked as germline (will not include variants QC) |
| `germline.vcf.gz.tbi` | Index file for germline small variant calls                           |
| `indel.vcf.gz`        | Raw indel calls in vcf format                                         |
| `indel.vcf.gz.tbi`    | Index for somatic indel calls                                         |
| `snv.vcf.gz`          | Raw SNV calls in vcf format                                           |
| `snv.vcf.gz.tbi`      | Index for SNV calls                                                   |
| `somatic.vcf.gz`      | SNV and indel calls marked as PASS and without a germline tag         |
| `somatic.vcf.gz.tbi`  | Index for somatic small variant calls                                 |

#### `severus`

```
├── severus
│   ├── all_SVs
│   │   ├── plots
│   │   │   ├── severus_{*}.html
│   │   ├── breakpoint_cluster_list.tsv
│   │   ├── breakpoint_clusters.tsv
│   │   ├── severus_all.vcf.gz
│   │   ├── severus_all.vcf.gz.tbi
│   ├── somatic_SVs
│   │   ├── plots
│   │   │   ├── severus_{*}.html
│   │   ├── breakpoint_cluster_list.tsv
│   │   ├── breakpoint_clusters.tsv
│   │   ├── severus_somatic.vcf.gz
│   │   ├── severus_somatic.vcf.gz.tbi
│   ├── breakpoints_double.csv
│   ├── read_ids.csv
│   ├── read_qual.txt
│   ├── severus.log
```

| File                                      | Description                                                                       |
| ----------------------------------------- | --------------------------------------------------------------------------------- |
| `all_SVs/plots/severus_{*}.html`          | html file containing a plot of connected breakpoints in a cluster                 |
| `all_SVs/breakpoint_cluster_list.tsv`     | tsv containing the breakpoints in all clustered events                            |
| `all_SVs/breakpoint_cluster.tsv`          | a tsv containing all clustered events                                             |
| `all_SVs/severus_all.vcf.gz`              | A vcf file containing all identified structural variants                          |
| `somatic_SVs/plots/severus_{*}.html`      | html file containing a plot of connected breakpoints in a cluster                 |
| `somatic_SVs/breakpoint_cluster_list.tsv` | tsv containing the breakpoints in somatic clustered events                        |
| `somatic_SVs/breakpoint_cluster.tsv`      | a tsv containing somatic clustered events                                         |
| `somatic_SVs/severus_somatic.vcf.gz`      | A vcf file containing identified somatic structural variants                      |
| `somatic_SVs/severus_somatic.vcf.gz.tbi`  | Index for identified somatic structural variants                                  |
| `breakpoints_double.csv`                  | csv file containing detailed information about identified breakpoints in bam file |
| `read_ids.csv`                            | a csv file containing read ids associated with each identified SV                 |
| `read_qual.txt`                           | file containing quality statistics about identified segements                     |
| `severus.log`                             | log file                                                                          |

#### `savana`

SAVANA structural variant and copy-number calling. Runs alongside Severus/ASCAT rather than replacing
either. Matched tumor/normal samples run `savana run` + `savana classify` + `savana cna` as separate
steps; tumor-only samples run the combined `savana to` command instead, producing the same file set
in one step. We strongly recommend matched tumor/normal mode for best performance -- SAVANA's own
docs note tumor-only calling is a fallback, best combined with population/panel-of-normals filtering.

```
├── savana
│   ├── sample.sv_breakpoints.vcf
│   ├── sample.sv_breakpoints.bedpe
│   ├── sample.sv_breakpoints_read_support.tsv
│   ├── sample.inserted_sequences.fa
│   ├── sample.classified.vcf
│   ├── sample.classified.somatic.vcf
│   ├── sample.classified.somatic.bedpe
│   ├── sample.classified.germline.vcf
│   ├── sample_segmented_absolute_copy_number.tsv
│   ├── sample_ranked_solutions.tsv
│   ├── sample_fitted_purity_ploidy.tsv
│   ├── sample_raw_read_counts.tsv
│   ├── sample_read_counts_<mode>_log2r_segmented.tsv
│   ├── sample_allele_counts_hetSNPs.bed
│   ├── 10kbp_bin_ref_all_sample.bed
```

| File                                          | Description                                                                     |
| ---------------------------------------------- | -------------------------------------------------------------------------------- |
| `sample.sv_breakpoints.vcf`                    | Raw (unclassified) SV breakpoints from `savana run`/`savana to`                  |
| `sample.sv_breakpoints.bedpe`                  | Raw SV breakpoints in BEDPE format                                               |
| `sample.sv_breakpoints_read_support.tsv`       | Supporting-read evidence per breakpoint                                          |
| `sample.inserted_sequences.fa`                 | Inserted sequences at breakpoints (insertion SVs)                                |
| `sample.classified.vcf`                        | All breakpoints after `savana classify` (somatic + germline)                     |
| `sample.classified.somatic.vcf`                | Classified somatic SV VCF -- fed into VEP for annotation                         |
| `sample.classified.somatic.bedpe`              | Classified somatic SVs in BEDPE format                                           |
| `sample.classified.germline.vcf`               | Classified germline SVs                                                          |
| `sample_segmented_absolute_copy_number.tsv`    | Segmented absolute copy-number calls from `savana cna`/`savana to`               |
| `sample_ranked_solutions.tsv`                  | Candidate purity/ploidy solutions, ranked                                        |
| `sample_fitted_purity_ploidy.tsv`              | Selected purity/ploidy fit                                                       |
| `sample_raw_read_counts.tsv`                   | Raw binned read counts used for CN segmentation                                  |
| `sample_read_counts_<mode>_log2r_segmented.tsv`| Segmented log2 read-count ratios per bin (`<mode>` is `mnorm` or `self`)          |
| `sample_allele_counts_hetSNPs.bed`             | Heterozygous-SNP allele counts (only when SNP/allele-frequency input is given)    |
| `10kbp_bin_ref_all_sample.bed`                 | Genome bins used for read-count binning                                          |

#### `deepvariant`

DeepVariant germline small variant calls. Present in all samples.

```
├── deepvariant
│   ├── sample.vcf.gz
│   ├── sample.vcf.gz.tbi
│   ├── sample.g.vcf.gz        # only when --generate_gvcf is true
│   ├── sample.g.vcf.gz.tbi    # only when --generate_gvcf is true
```

| File                  | Description                                                                     |
| --------------------- | ------------------------------------------------------------------------------- |
| `sample.vcf.gz`       | DeepVariant germline SNV and indel calls in VCF format                          |
| `sample.vcf.gz.tbi`   | Index for DeepVariant germline calls                                            |
| `sample.g.vcf.gz`     | DeepVariant gVCF file with calls at all positions (only with `--generate_gvcf`) |
| `sample.g.vcf.gz.tbi` | Index for DeepVariant gVCF (only with `--generate_gvcf`)                        |

#### `deepsomatic`

DeepSomatic somatic small variant calls. Present in all samples.

```
├── deepsomatic
│   ├── sample.vcf.gz
│   ├── sample.vcf.gz.tbi
│   ├── sample.g.vcf.gz        # only when --generate_gvcf is true
│   ├── sample.g.vcf.gz.tbi    # only when --generate_gvcf is true
```

| File                  | Description                                                                     |
| --------------------- | ------------------------------------------------------------------------------- |
| `sample.vcf.gz`       | DeepSomatic somatic SNV and indel calls in VCF format                           |
| `sample.vcf.gz.tbi`   | Index for DeepSomatic somatic calls                                             |
| `sample.g.vcf.gz`     | DeepSomatic gVCF file with calls at all positions (only with `--generate_gvcf`) |
| `sample.g.vcf.gz.tbi` | Index for DeepSomatic gVCF (only with `--generate_gvcf`)                        |

#### `phased`

Phased variant calls produced by Longphase. Present in all samples.

```
├── phased
│   ├── germline_smallvariants.vcf.gz
│   ├── germline_smallvariants.vcf.gz.tbi
│   ├── somatic_smallvariants.vcf.gz
│   ├── somatic_smallvariants.vcf.gz.tbi
```

| File                                | Description                                                      |
| ----------------------------------- | ---------------------------------------------------------------- |
| `germline_smallvariants.vcf.gz`     | Longphase-phased germline SNV/indel VCF with haplotype (PS) tags |
| `germline_smallvariants.vcf.gz.tbi` | Index for the phased germline VCF                                |
| `somatic_smallvariants.vcf.gz`      | Longphase-phased somatic SNV/indel VCF with haplotype (PS) tags  |
| `somatic_smallvariants.vcf.gz.tbi`  | Index for the phased somatic VCF                                 |

</details>

### `vep`

<details markdown="1">
<summary>Output files</summary>

```
├── vep
│   ├── germline
│   │   ├── sample_GERMLINE_VEP.vcf.gz
│   │   ├── sample_GERMLINE_VEP_summary.html
│   │   ├── sample_GERMLINE_VEP.vcf.gz.tbi
│   ├── somatic
│   │   ├── sample_SOMATIC_VEP.vcf.gz
│   │   ├── sample_SOMATIC_VEP_summary.html
│   │   ├── sample_SOMATIC_VEP.vcf.gz.tbi
│   ├── SVs
│   │   ├── sample_SV_VEP.vcf.gz
│   │   ├── sample_SV_VEP_summary.html
│   │   ├── sample_SV_VEP.vcf.gz.tbi
```

| File                                        | Description                                                             |
| ------------------------------------------- | ----------------------------------------------------------------------- |
| `germline/sample_GERMLINE_VEP.vcf.gz`       | Annotated germline indel and SNV vcf file                               |
| `germline/sample_GERMLINE_VEP_summary.html` | Visual summary of germline indel and SNV annotations in html format     |
| `germline/sample_GERMLINE_VEP.vcf.gz.tbi`   | Annotated germline indel and SNV vcf index file                         |
| `somatic/sample_SOMATIC_VEP.vcf.gz`         | Annotated somatic indel and SNV vcf file                                |
| `somatic/sample_SOMATIC_VEP_summary.html`   | Visual summary of somatic indel and SNV annotations in html format      |
| `somatic/sample_SOMATIC_VEP.vcf.gz.tbi`     | Annotated somatic indel and SNV vcf index file                          |
| `SVs/sample_SV_VEP.vcf.gz`                  | Annotated somatic structural variant vcf file                           |
| `SVs/sample_SV_VEP_summary.html`            | Visual summary of somatic structural variant annotations in html format |
| `SVs/sample_SV_VEP.vcf.gz.tbi`              | Annotated somatic structural variant vcf index file                     |

</details>

### `wakhan`

<details markdown="1">
<summary>Output files</summary>

```
├── wakhan
│   ├── {ploidy}_{purity}_{confidence}
│   │   ├── bed_output
│   │   │   ├── genes_copynumber_states.bed
│   │   │   ├── loh_regions.bed
│   │   │   ├── sample_{ploidy}_{purity}_{confidence}_HP_1.bed
│   │   │   ├── sample_{ploidy}_{purity}_{confidence}_HP_2.bed
│   │   ├── variation_plots
│   │   │   ├── chr{1-22,X,Y}_cn.html
│   │   │   ├── chr{1-22,X,Y}_cn.pdf
│   │   │   ├── CN_VARIATION_INDEX.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genes_genome.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genes_genome.pdf
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_breakpoints.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_breakpoints.pdf
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.html
│   │   ├── sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.pdf
│   ├── coverage_data
│   │   ├── {0-23}_SNPS.csv
│   │   ├── coverage_ps.csv
│   │   ├── phase_corrected_coverage.csv
│   │   ├── pileup_SNPs.csv
│   ├── coverage_plots
│   │   ├── chr{1-22,X,Y}_cov.html
│   │   ├── chr{1-22,X,Y}.pdf
│   │   ├── COVERAGE_INDEX.html
│   ├── phasing_output
│   │   ├── chr{1-22,X,Y}_phase_correction_0.html
│   │   ├── chr{1-22,X,Y}_phase_correction_1.html
│   │   ├── chr{1-22,X,Y}_without_phase_correction.html
│   │   ├── chr{1-22,X,Y}.pdf
│   │   ├── sample.rephased.vcf.gz
│   │   ├── sample.rephased.vcf.gz.tbi
│   ├── snps_loh_plots
│   │   ├── chr{1-22,X,Y}_snps_loh.html
│   ├── sample_heatmap_ploidy_purity.html
│   ├── sample_heatmap_ploidy_purity.html.pdf
│   ├── sample_optimized_peak.html
│   ├── solutions_ranks.tsv

```

| File                                                                                                   | Description                                                                                        |
| ------------------------------------------------------------------------------------------------------ | -------------------------------------------------------------------------------------------------- |
| `{ploidy}_{purity}_{confidence}/bed_output/genes_copynumber_states.bed`                                | bed file containing allele specific copy number values with coverage information                   |
| `{ploidy}_{purity}_{confidence}/bed_output/loh_regions.bed`                                            | bed file containing positions of loss of heterozygosity regions                                    |
| `{ploidy}_{purity}_{confidence}/bed_output/sample_{ploidy}_{purity}_{confidence}_HP_1.bed`             | bed file containing copy number states, coverage, and SV breakpoints for haplotype 1               |
| `{ploidy}_{purity}_{confidence}/bed_output/sample_{ploidy}_{purity}_{confidence}_HP_2.bed`             | bed file containing copy number states, coverage, and SV breakpoints for haplotype 2               |
| `{ploidy}_{purity}_{confidence}/variation_plots/chr{1-22,X,Y}_cn.html`                                 | html based plotly plot of copy number and coverage for individual chromosomes                      |
| `{ploidy}_{purity}_{confidence}/variation_plots/chr{1-22,X,Y}_cn.pdf`                                  | pdf based plotly plot of copy number and coverage for individual chromosomes                       |
| `{ploidy}_{purity}_{confidence}/variation_plots/CN_VARIATION_INDEX.html`                               | unclear html plot                                                                                  |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genes_genome.html`               | html plots of copy number variations in highlighted genes                                          |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genes_genome.pdf`                | pdf plots of copy number variations in highlighted genes                                           |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.html` | genome-wide html copy number plots with coverage information on same axis                          |
| `{ploidy}_{purity}_{confidence}/sample_{purity}_{ploidy}_{confidence}_genome_copynumbers_details.pdf`  | genome-wide pdf copy number plots with coverage information on same axis                           |
| `coverage_data/{0-23}_SNP.csv`                                                                         | CSV of coverage data per chromosome                                                                |
| `coverage_data/coverage_ps.csv`                                                                        | CSV of overall haplotype specific coverage data                                                    |
| `coverage_data/coverage.csv`                                                                           | CSV of overall coverage data                                                                       |
| `coverage_data/phase_corrected_coverage.csv`                                                           | CSV of overall phase-corrected coverage data                                                       |
| `coverage_data/pileup_SNPs.csv`                                                                        | CSV of SNP pileup data                                                                             |
| `coverage_plots/chr{1-22,X,Y}_cov.html`                                                                | chromosome specific html coverage plots                                                            |
| `coverage_plots/chr{1-22,X,Y}_cov.pdf`                                                                 | chromosome specific pdf coverage plots                                                             |
| `coverage_plots/COVERAGE_INDEX.html`                                                                   | unclear html plot                                                                                  |
| `phasing_output/chr{1-23,X,Y}_phase_correction_0.html`                                                 | Phase-switch error correction plot per chromosome                                                  |
| `phasing_output/chr{1-23,X,Y}_phase_correction_1.html`                                                 | Phase-switch error correction plot per chromosome                                                  |
| `phasing_output/chr{1-22,X,Y}_without_phase_correction.html`                                           | Phase-switch error without phase correction plot per chromosome                                    |
| `phasing_output/chr{1-22,X,Y}.pdf`                                                                     | Phase-switch error correction plot                                                                 |
| `phasing_output/sample_rephased.vcf.gz`                                                                | phase corrected SNP vcf file                                                                       |
| `phasing_output/sample_rephased.vcf.gz.tbi`                                                            | phase corrected SNP vcf index file                                                                 |
| `snps_loh_plots/chr{1-22,X,Y}_snps_loh.html`                                                           | interactive HTML plots of SNP allele frequencies and loss of heterozygosity regions per chromosome |
| `sample_heatmap_ploidy_purity.html`                                                                    | heatmap html plot of purity ploidy fit                                                             |
| `sample_heatmap_ploidy_purity.html.pdf`                                                                | heatmap pdf plot of purity ploidy fit                                                              |
| `sample_optimized_peak.html`                                                                           | optimization peak plot                                                                             |
| `solutions_ranks.tsv`                                                                                  | rank of potential purity ploidy solutions                                                          |

</details>

### `multiqc`

<details markdown="1">
<summary>Output files</summary>

```
├── multiqc
│   ├── multiqc_data
│   │   ├── BETA-multiqc.parquet
│   │   ├── llms-full.txt
│   │   ├── mosdepth_cov_dist.txt
│   │   ├── mosdepth_cumcov_dist.txt
│   │   ├── mosdepth_perchrom.txt
│   │   ├── mosdepth-coverage-per-contig-multi.txt
│   │   ├── mosdepth-cumcoverage-dist-id.txt
│   │   ├── mosdepth-xy-coverage-plot.txt
│   │   ├── multiqc_citations
│   │   ├── multiqc_data.json
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc_software_versions.txt
│   │   ├── multiqc_sources.txt
│   │   ├── multiqc.log
│   ├── multiqc_plots
│   │   ├── pdf
│   │   │   ├── mosdepth-coverage-per-contig-multi-cnt.pdf
│   │   │   ├── mosdepth-coverage-per-contig-multi-log.pdf
│   │   │   ├── mosdepth-cumcoverage-dist-id.pdf
│   │   │   ├── mosdepth-xy-coverage-plot-cnt.pdf
│   │   │   ├── mosdepth-xy-coverage-plot-pct.pdf
│   │   ├── png
│   │   │   ├── mosdepth-coverage-per-contig-multi-cnt.png
│   │   │   ├── mosdepth-coverage-per-contig-multi-log.png
│   │   │   ├── mosdepth-cumcoverage-dist-id.png
│   │   │   ├── mosdepth-xy-coverage-plot-cnt.png
│   │   │   ├── mosdepth-xy-coverage-plot-pct.png
│   │   ├── svg
│   │   │   ├── mosdepth-coverage-per-contig-multi-cnt.svg
│   │   │   ├── mosdepth-coverage-per-contig-multi-log.svg
│   │   │   ├── mosdepth-cumcoverage-dist-id.svg
│   │   │   ├── mosdepth-xy-coverage-plot-cnt.svg
│   │   │   ├── mosdepth-xy-coverage-plot-pct.svg
│   ├── multiqc_report.html

```

| File                                                           | Description                                                  |
| -------------------------------------------------------------- | ------------------------------------------------------------ |
| `multiqc_data/BETA-multiqc.parquet`                            | Multiqc data in Apache Parquet format (BETA)                 |
| `multiqc_data/llms-full.txt`                                   | Prompt for large-language-model summary                      |
| `multiqc_data/mosdepth_cov_dist.txt`                           | text file of coverage distribution                           |
| `multiqc_data/mosdepth_cumcov_dist.txt`                        | text file of cummulative coverage distribution               |
| `multiqc_data/mosdepth_perchrom.txt`                           | text file of coverage per chromosome                         |
| `multiqc_data/mosdepth-coverage-per-contig-multi.txt`          | text file of coverage per contig                             |
| `multiqc_data/mosdepth-cumcoverage-dist-id.txt`                | unclear text file                                            |
| `multiqc_data/mosdepth-xy-coverage-plot.txt`                   | summary of chr X, chr Y coverage                             |
| `multiqc_data/multiqc_citations`                               | citations for multiqc                                        |
| `multiqc_data/multiqc_data.json`                               | json file containing multiqc data output                     |
| `multiqc_data/multiqc_general_stats.txt`                       | summary statistics for the samples                           |
| `multiqc_data/multiqc_software_versions.txt`                   | software versions for tools used in multiqc                  |
| `multiqc_data/multiqc_sources.txt`                             | software sources for tools used in multiqc                   |
| `multiqc_data/multiqc.log`                                     | log file for multiqc                                         |
| `multiqc_plots/pdf/mosdepth-coverage-per-contig-multi-cnt.pdf` | pdf format plot of coverage per contig                       |
| `multiqc_plots/pdf/mosdepth-coverage-per-contig-multi-log.pdf` | pdf format plot of coverage per contig on a logarithmic plot |
| `multiqc_plots/pdf/mosdepth-cumcoverage-dist-id.pdf`           | pdf format plot of distribution of cumulative coverage       |
| `multiqc_plots/pdf/mosdepth-xy-coverage-plot-cnt.pdf`          | pdf format plot of chr X and chr Y coverage by count         |
| `multiqc_plots/pdf/mosdepth-xy-coverage-plot-pct.pdf`          | pdf format plot of chr X and chr Y coverage by percentage    |
| `multiqc_plots/png/mosdepth-coverage-per-contig-multi-cnt.png` | png format plot of coverage per contig                       |
| `multiqc_plots/png/mosdepth-coverage-per-contig-multi-log.png` | png format plot of coverage per contig on a logarithmic plot |
| `multiqc_plots/png/mosdepth-cumcoverage-dist-id.png`           | png format plot of distribution of cumulative coverage       |
| `multiqc_plots/png/mosdepth-xy-coverage-plot-cnt.png`          | png format plot of chr X and chr Y coverage by count         |
| `multiqc_plots/png/mosdepth-xy-coverage-plot-pct.png`          | png format plot of chr X and chr Y coverage by percentage    |
| `multiqc_plots/svg/mosdepth-coverage-per-contig-multi-cnt.svg` | svg format plot of coverage per contig                       |
| `multiqc_plots/svg/mosdepth-coverage-per-contig-multi-log.svg` | svg format plot of coverage per contig on a logarithmic plot |
| `multiqc_plots/svg/mosdepth-cumcoverage-dist-id.svg`           | svg format plot of distribution of cumulative coverage       |
| `multiqc_plots/svg/mosdepth-xy-coverage-plot-cnt.svg`          | svg format plot of chr X and chr Y coverage by count         |
| `multiqc_plots/svg/mosdepth-xy-coverage-plot-pct.svg`          | svg format plot of chr X and chr Y coverage by percentage    |
| `multiqc_plots/multiqc_report.html`                            | svg format plot of chr X and chr Y coverage by percentage    |

### `pipeline_info`

<details markdown="1">
<summary>Output files</summary>

```
├── pipeline_info
│   ├── execution_timeline_{DATE}.html
│   ├── execution_trace_{DATE}.txt
│   ├── final_sample_disk_usage.tsv
│   ├── lrsomatic_software_mqc_versions.yml
│   ├── params_{DATE}.json
│   ├── pipeline_dag_{DATE}.html
│   ├── raw_task_disk_usage.tsv
```

| File                                  | Description                                                                                 |
| ------------------------------------- | ------------------------------------------------------------------------------------------- |
| `execution_timeline_{DATE}.html`      | a graphical summary of the timing of each module's task over the course of the pipeline run |
| `execution_trace_{DATE}.txt`          | detailed per-task resource usage log (CPU, memory, wall time)                               |
| `final_sample_disk_usage.tsv`         | summary of disk usage per sample at pipeline completion                                     |
| `lrsomatic_software_mqc_versions.yml` | summary of the versions of each tool used by the pipeline                                   |
| `params_{DATE}.json`                  | summary of the parameters used in the pipeline run                                          |
| `pipeline_dag_{DATE}.html`            | flow chart summarizing the pipeline structure                                               |
| `raw_task_disk_usage.tsv`             | per-task disk usage across all pipeline tasks                                               |

</details>
