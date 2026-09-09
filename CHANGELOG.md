# IntGenomicsLab/lrsomatic: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev

### `Added`

- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added `LRSOMATICREPORT` as the final pipeline step: a self-contained per-sample HTML report covering small variants, structural variants, copy number and QC. Skip it with `--skip_report`; choose the gene panel selected on load with `--report_gene_panel` (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Vendored the [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) v1.3.2 tool source at `assets/lrsomatic_report`, so `nextflow run IntGenomicsLab/lrsomatic` ships it without a submodule checkout (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Added a `solution_dirs` output to the WAKHAN module so its per-solution copy-number plots can be staged downstream (@ljwharbers).

### `Changed`

- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - Re-synced the vendored [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) to v1.3.0, which adds tickbox dropdown filters on the categorical columns of both variant tables and turns the report's gene panel selector into checkboxes (@ljwharbers).
- [#176](https://github.com/IntGenomicsLab/lrsomatic/pull/176) - Re-synced the vendored [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) to v1.3.2: facet dropdown counts follow the active filters, opening a facet menu no longer resets the table's horizontal scroll, a flatter clinical theme, and inline code comments trimmed to one line (@ljwharbers).
- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - `--report_gene_panel` now takes a comma-separated list, so several panels can be applied at once: a variant or SV is kept if it hits any of them. Panel values are also validated at launch instead of failing inside the report task (@ljwharbers).
- [#184](https://github.com/IntGenomicsLab/lrsomatic/pull/184) - Replaced the CHM13 Severus panel of normals with the merged 1000 Genomes + ASAP panel (@AmberVerhasselt).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Updated to the nf-core template for [nf-core/tools 4.1.0](https://github.com/nf-core/tools/releases/tag/4.1.0). Adds the `check-added-large-files`, `check-merge-conflict` and `block-pipeline-outdir` pre-commit hooks, a `process_low_memory` resource label, and a split-out `pr-comment.yml` workflow. `conf/igenomes.config` was converted to the template's strict-syntax `params.genomes` map literal, and the MultiQC module was bumped to 1.35 (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Raised the minimum Nextflow version to `25.10.4`, following the template, and bumped the nf-test CI matrix to match (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Filled in the `manifest.contributors` `contribution` fields: @ljwharbers and @robert-a-forsyth as author and maintainer, all other contributors as author. Previously empty, which left the RO-Crate metadata without any authors (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Added @AmberVerhasselt to `manifest.contributors` as a contributor and to the README credits, and expanded @ljwharbers' affiliation to match the other core contributors. RO-Crate metadata updated to match (@ljwharbers).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Set @laulambr's affiliation to the same three-part VIB/KU Leuven string as the other core contributors (@ljwharbers).

### `Fixed`

- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - NanoPlot pre- and post-alignment statistics now reach MultiQC. `ch_nanoplot_pre_txt` and `ch_nanoplot_post_txt` were declared empty and mixed into the MultiQC inputs, but never assigned from `NANOPLOT_PRE.out.txt` / `NANOPLOT_POST.out.txt`, so the NanoStat section was silently missing from every report (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - `NANOPLOT_PRE` now includes the replicate in its output prefix, so samples with more than one replicate no longer collapse into a single MultiQC sample (@ljwharbers).
- [#181](https://github.com/IntGenomicsLab/lrsomatic/pull/181) - samtools stats/flagstat/idxstats and mosdepth outputs are now prefixed `<sample>_<type>` instead of `<sample>`. The tumor and normal files of a matched pair shared a basename, so MultiQC logged `Duplicate sample name found! Overwriting` and reported only one of the two; the general statistics table now has one row per BAM (`<sample>_tumor`, `<sample>_normal`), which the post-alignment NanoStat row joins (@ljwharbers).
- [#186](https://github.com/IntGenomicsLab/lrsomatic/pull/186) - Stopped snapshotting the md5 of sample4's merged tumour BAM and its index in the `clair_only` nf-test: `samtools merge` gives the colliding `@PG` IDs of the two replicates a random hex suffix, so neither digest is reproducible. The alignment records are, and are now asserted with `bam().getReadsMD5()` instead (@ljwharbers).
- [#182](https://github.com/IntGenomicsLab/lrsomatic/pull/182) - Added `--vcf` to the default `vep_args` so VEP writes VCF output rather than its default tab-delimited format (@AmberVerhasselt).
- [#183](https://github.com/IntGenomicsLab/lrsomatic/pull/183) - Corrected the `github` URL for Laurens Lambrechts in `manifest.contributors`, which was a copy of @MariosEft97's, to @laulambr. RO-Crate metadata updated to match (@ljwharbers).

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
