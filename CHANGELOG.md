# IntGenomicsLab/lrsomatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev

### `Added`

- [#189](https://github.com/IntGenomicsLab/lrsomatic/pull/189) - Added SAVANA structural variant and copy-number calling, running alongside Severus/ASCAT (@yannvrb).

### `Changed`

- [#184](https://github.com/IntGenomicsLab/lrsomatic/pull/184) - Replaced the CHM13 Severus panel of normals with the merged 1000 Genomes + ASAP panel (@AmberVerhasselt).

### `Fixed`

- [#182](https://github.com/IntGenomicsLab/lrsomatic/pull/182) - Added `--vcf` to the default `vep_args` so VEP writes VCF output rather than its default tab-delimited format (@AmberVerhasselt).

## v1.1.0 - [2026-04-28]

### `Added`

- [#117](https://github.com/IntGenomicsLab/lrsomatic/pull/117) - Added ASCAT PDF plots to output (@robert-a-forsyth).
- [#126](https://github.com/IntGenomicsLab/lrsomatic/pull/126) - Added ASCAT raw segments txt files to output (@AmberVerhasselt).
- [#135](https://github.com/IntGenomicsLab/lrsomatic/pull/135) - Added `skip_m6a` parameter to allow skipping m6A base modification steps (@robert-a-forsyth).
- [#141](https://github.com/IntGenomicsLab/lrsomatic/pull/141) - Added output of phased variants in separate VCF files for improved downstream analysis (@ljwharbers).
- [#143](https://github.com/IntGenomicsLab/lrsomatic/pull/143) - Added Severus `min_support` parameter and `skip_fibernormal` option (@ljwharbers).
- [#145](https://github.com/IntGenomicsLab/lrsomatic/pull/145) - Integrated MultiQC and nanoplot for comprehensive QC reporting with long-read sequencing metrics (@ljwharbers).
- [#147](https://github.com/IntGenomicsLab/lrsomatic/pull/147) - Implemented whatshap_stats module to generate phase block statistics and phasing quality metrics (@ljwharbers).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Added DeepVariant and DeepSomatic modules for germline and somatic variant calling from long-read sequencing data (@robert-a-forsyth).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Added GPU support for Clair3, DeepVariant, and fibertools (@robert-a-forsyth).
- [#150](https://github.com/IntGenomicsLab/lrsomatic/pull/150) - Added Claude GitHub Actions workflows for automated code review and PR assistance (@ljwharbers).
- [#152](https://github.com/IntGenomicsLab/lrsomatic/pull/152) - Integrated modkit module for long-read base modification detection and analysis (@robert-a-forsyth).
- [#165](https://github.com/IntGenomicsLab/lrsomatic/pull/165) - Added bcftools/view and samtools/merge modules; added extended test suites for union, consensus, clair-only, and deep-only caller modes (@robert-a-forsyth).

### `Changed`

- [#123](https://github.com/IntGenomicsLab/lrsomatic/pull/123) - Updated channel structure (@robert-a-forsyth).
- [#137](https://github.com/IntGenomicsLab/lrsomatic/pull/137) - Bulk module versions update. Fixed some issues with Wakhan (@ljwharbers).
- [#138](https://github.com/IntGenomicsLab/lrsomatic/pull/138) - Perform QC before merging replicates (@robert-a-forsyth).
- [#140](https://github.com/IntGenomicsLab/lrsomatic/pull/140) - Improved documentation with additional pipeline usage examples and configuration guidance (@ljwharbers).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Refactored variant calling workflow to support both DeepVariant and existing callers with improved configuration handling (@robert-a-forsyth).
- [#152](https://github.com/IntGenomicsLab/lrsomatic/pull/152) - Updated container versions and dependencies for modkit and related tools (@robert-a-forsyth).
- [#157](https://github.com/IntGenomicsLab/lrsomatic/pull/157) - Added ASAP Panel of Normals citation to CITATIONS.md (@ljwharbers).
- [#160](https://github.com/IntGenomicsLab/lrsomatic/pull/160) - DeepVariant/DeepSomatic optimization and Panel-of-Normals handling improvements; updated docs; removed Claude workflow files (@robert-a-forsyth).
- [#163](https://github.com/IntGenomicsLab/lrsomatic/pull/163) - Added LongPhase supplementary alignment tag to extended args (@AmberVerhasselt).
- [#164](https://github.com/IntGenomicsLab/lrsomatic/pull/164) - Updated nf-core template components (@robert-a-forsyth).
- [#166](https://github.com/IntGenomicsLab/lrsomatic/pull/166) - Updated Wakhan to v0.4.3 using BioContainers distribution (@robert-a-forsyth).

### `Fixed`

- [#116](https://github.com/IntGenomicsLab/lrsomatic/pull/116) - Corrected ASCAT GC and RT bias correction (@AmberVerhasselt).
- [#118](https://github.com/IntGenomicsLab/lrsomatic/pull/118) - Updated nf-core template components to align with latest pipeline standards (@ljwharbers).
- [#130](https://github.com/IntGenomicsLab/lrsomatic/pull/130) - Fixed path pattern for rephased VCF files (@Tim-Yu).
- [#137](https://github.com/IntGenomicsLab/lrsomatic/pull/137) - Resolved Nextflow strict syntax compliance issues for compatibility with latest Nextflow versions (@ljwharbers).
- [#149](https://github.com/IntGenomicsLab/lrsomatic/pull/149) - Corrected bcftools and vcfsplit operations for accurate variant filtering and merging (@robert-a-forsyth).
- [#165](https://github.com/IntGenomicsLab/lrsomatic/pull/165) - Fixed consensus variant calling workflow issues (@robert-a-forsyth).
- [#168](https://github.com/IntGenomicsLab/lrsomatic/pull/168) - Fixed default values for `germline_var_keep`, `somatic_var_keep`, `prioritize_caller_germline`, and `prioritize_caller_somatic` parameters to default to `clair` (@robert-a-forsyth).

## v1.0.0 - [28 Nov 2025]

Initial release of IntGenomicsLab/lrsomatic, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- ClairS-TO Module
- cramino module
- modkit pileup module
- mosdepth module
- minimap2/index module
- minimap2/align module
- pigz module
- samtools/cat module
- samtools faidx module
- bam_stats_samtools subworkflow
- mosdepth added to workflow
- add longphase/tag and longphase/phase modules

### `Fixed`

- New channel structure
- No longer possible to have duplicated naming after samtools cat
- restructured minimap2

### `Dependencies`

### `Deprecated`
