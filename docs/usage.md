# IntGenomicsLab/lrsomatic: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file of the following form:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber
sample1,tumour.bam,normal.bam,ont,female,n
sample2,tumour.bam,,ont,female,y
sample3,tumour.bam,,pb,male,n
sample4,tumour.bam,normal.bam,pb,male,y
```

lrsomatic extracts information from the bam header files to decide which models to use for Clair3, ClairS, or ClairS-TO. However, this can optionally be specified manually. You can do this for one or many samples, if the field is left blank, the pipeline will default to extracting this information. You can specify this by creating your csv in the following form:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber,clair3_model,clairSTO_model,clairS_model
sample1,tumour.bam,normal.bam,ont,female,n
sample2,tumour.bam,,ont,female,y
sample3,tumour.bam,normal.bam,pb,male,n,r1041_e82_400bps_sup_v420,,ont_r10_dorado_sup_5khz_ssrs
sample4,tumour.bam,normal.bam,pb,male,y
```

Use the `input` parameter to specify the location to this input csv.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber
sample1,tumour1.bam,normal.bam,ont,female,n
sample1,tumour2.bam,,ont,female,y
sample1,tumour3.bam,,pb,male,n
```

### Full Description of Samplesheet Columns

| Column           | Description                                                                                                                                                                            |
| ---------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`         | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `bam_tumor`      | Full path to BAM file for the tumor. File must end in `.bam`.                                                                                                                          |
| `bam_normal`     | Full path to BAM file for the tumor. File must end in `.bam`.                                                                                                                          |
| `platform`       | A string specifying the platform used for sequencing, can be either `pb` for PacBio sequencing data or `ont` for Oxford Nanopore sequencing data                                       |
| `sex`            | A string specifying the biological sex of the sample, can either be `m` or `f`                                                                                                         |
| `fiber`          | A string specifying if the sample has been subjected to Fiber-seq. Can either be `y` or `n`                                                                                            |
| `clair3_model`   | A string describing which model is to be used for Clair3's small variant calling (_optional_)                                                                                          |
| `clairSTO_model` | A string describing which model is to be used for ClairS-TO's small variant calling (_optional_)                                                                                       |
| `clairS_model`   | A string describing which model is to be used for ClairS's small variant calling (_optional_)                                                                                          |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run IntGenomicsLab/lrsomatic --input ./samplesheet.csv --outdir ./results --genome GRCh38 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.

```

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run IntGenomicsLab/lrsomatic -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

## CHM13 Support

Our pipeline fully supports CHM13 and most reference and annotation files are automatically downloaded when specifying `--genome CHM13`.

However, VEP will need a bit of additional setup. The VEP cache for CHM13 needs to be manually downloaded. This can be done using the following code. Feel free to change any of the paths, ensuring that the correct path is pointed to in the pipeline parameters.

Download CHM13 Cache:

```bash
cd $HOME/.vep
curl -O https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/variation/2022_10/indexed_vep_cache/Homo_sapiens-GCA_009914755.4-2022_10.tar.gz
tar xzf Homo_sapiens-GCA_009914755.4-2022_10.tar.gz
```

Then you can run the pipeline as follows:

```bash
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv \
  --outdir ./results \
  --genome CHM13 \
  --vep_cache $HOME/.vep \
  --vep_cache_version 107 \
  -profile docker
```

If you want to run with a CHM13 reference without using `--genome CHM13` (for example, via a custom FASTA or configuration), you must also specify `--vep_genome T2T-CHM13v2.0` and `--vep_species homo_sapiens_gca009914755v4`.

This cache carries gene and transcript models only, so every pathogenicity score on CHM13 has to come from a plugin. With `--genome CHM13` the pipeline fetches those itself; note that only the two predictors keyed on protein rather than genome coordinates can reach the assembly at all. See [VEP plugins](#vep-plugins) for which those are and why.

Setting `--vep_genome T2T-CHM13v2.0` by hand, without `--genome CHM13`, resolves no plugin defaults — supply the `--vep_*` paths yourself, or accept VEP without plugin annotation.
For mutational signatures, `--genome CHM13` selects the `CHM13-T2T` SigProfilerMatrixGenerator genome. Its payload is not on the AlexandrovLab FTP yet, so `--download_sigprofiler_genome` fetches it from the IntGenomicsLab Globus collection (`--sigprofiler_genome_url`); see [Mutational Signature Options](#mutational-signature-options).

For structural variants, the CHM13 panel of normals is a merged panel combining the 1000 Genomes CHM13 panel shipped with SEVERUS and the ASAP cohort, with median confidence intervals per breakpoint. The pipeline exposes it as `--pon_file` and hands it to SEVERUS via that tool's own `--PON` flag; it is downloaded automatically with `--genome CHM13`. GRCh38 continues to use the 1000 Genomes panel shipped with SEVERUS.

### Pipeline options

| Parameter        | Description                                                                                                                                                                  |
| ---------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-input`         | Full file path to input samplesheet, must be in `.csv` format and conform to specifications noted above                                                                      |
| `--genome`       | Specified genome assembly, support is given for `GRCh38` and `CHM13`                                                                                                         |
| `--normal_fiber` | A boolean which skips fiber-seq processing on normal files (on those which have fiber-seq for the tumor). Default = `true` (_does not skip fiber-seq processing for normal_) |

#### Skipping options:

| Parameter              | Description                                                                                                                                                                          |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--skip_qc`            | A boolean to skip all QC steps, including `mosdepth`, `samtools`,`fibertools`, `cramino`. Default = `false`                                                                          |
| `--skip_fiber`         | A boolean to skip all `fibertools` related modules. Default = `false`                                                                                                                |
| `--skip_cramino`       | A boolean to skip `cramino`. Default = `false`                                                                                                                                       |
| `--skip_mosdepth`      | A boolean to skip `mosdepth`. Default = `false`                                                                                                                                      |
| `--skip_ascat`         | A boolean to skip `ascat`. Default = `false`                                                                                                                                         |
| `--skip_bamstats`      | A boolean to skip `bamstats`. Default = `false`                                                                                                                                      |
| `--skip_wakhan`        | A boolean to skip `wakhan`. Default = `false`                                                                                                                                        |
| `--skip_vep`           | A boolean to skip `vep`. Default = `false`                                                                                                                                           |
| `--skip_m6a`           | A boolean to skip `fibertools_m6a`, used if you have m6a calls but would still like nucleosome positions for PacBio data (ONT data is required to have m6a calls). Default = `false` |
| `--skip_nanoplot`      | A boolean to skip NanoPlot QC on aligned and unaligned BAM files. Default = `false`                                                                                                  |
| `--skip_normalfiber`   | A boolean to skip fibertools processing for the normal sample. Default = `false`                                                                                                     |
| `--skip_modcall`       | A boolean to skip modkit methylation calling. Default = `false`                                                                                                                      |
| `--skip_modkit`        | A boolean to skip the modkit pileup step. Default = `false`                                                                                                                          |
| `--skip_whatshapstats` | A boolean to skip WhatsHap phasing statistics. Default = `false`                                                                                                                     |
| `--skip_signatures`    | A boolean to skip mutational signature analysis (SigProfilerMatrixGenerator + SigProfilerAssignment). Default = `false`                                                              |
| `--skip_report`        | A boolean to skip the final per-sample HTML report. Default = `false`                                                                                                                |

#### LONGPHASE options:

| Parameter                       | Description                                                                        |
| ------------------------------- | ---------------------------------------------------------------------------------- |
| `--longphase_tag_supplementary` | Include supplementary alignments in Longphase haplotype tagging. Default = `false` |

#### VEP options:

| Parameter              | Description                                                                                                                                      |
| ---------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--vep_cache`          | Full path to a vep cache. If left blank, this will default to pulling from this [Annotation Cache Storage](https://annotation-cache.github.io/). |
| `--vep_cache_version`  | Integer specifying version of vep cache. Default = `113`                                                                                         |
| `--vep_args`           | A string specifying arguments to vep. Default = `"--everything --filter_common --per_gene --total_length --offline --format vcf --vcf"`          |
| `--vep_custom`         | A full path to a vcf file containing custom variants for annotation. Must be bgzipped and have `.vcf.gz` format. Default = `null`                |
| `--vep_custom_tbi`     | A full path to a index file for cutom vcf for vep. Default = `null`                                                                              |
| `--download_vep_cache` | A boolean to automatically download the VEP cache if not found locally. Default = `false`                                                        |

#### VEP plugin options:

The plugin data is **on by default** on `--genome GRCh38` and `--genome CHM13`: each parameter
below falls back to a per-assembly default that the pipeline downloads and, where the release needs
it, reshapes. Setting one overrides that default; `--skip_vep_plugins` turns the whole set off. See
[VEP plugins](#vep-plugins) for sizes, licence terms and which assembly each applies to.

| Parameter                    | Description                                                                                                              |
| ---------------------------- | ------------------------------------------------------------------------------------------------------------------------ |
| `--skip_vep_plugins`         | Annotate with VEP alone, fetching no plugin data. Default = `false`                                                      |
| `--vep_alphamissense`        | AlphaMissense GRCh38 score file, for the `AlphaMissense` plugin. GRCh38 only                                             |
| `--vep_alphamissense_tbi`    | Index for `--vep_alphamissense`. Required whenever `--vep_alphamissense` is set                                          |
| `--vep_alphamissense_aa`     | AlphaMissense protein-space release, or a table already built from it, for the `AlphaMissenseProtein` plugin. CHM13 only |
| `--vep_alphamissense_aa_tbi` | Index for `--vep_alphamissense_aa`. Required whenever `--vep_alphamissense_aa` is set                                    |
| `--vep_polyphen_sift_db`     | Ensembl pangenome PolyPhen/SIFT SQLite database, for the `PolyPhen_SIFT` plugin. Needed on CHM13 only                    |
| `--vep_clinvar`              | ClinVar VCF, added as a VEP `--custom` annotation                                                                        |
| `--vep_clinvar_tbi`          | Index for `--vep_clinvar`. Required whenever `--vep_clinvar` is set                                                      |
| `--vep_clinvar_fields`       | Comma-separated ClinVar INFO fields to carry through. Default = `"CLNSIG,CLNREVSTAT,CLNDN"`                              |
| `--vep_cadd_snv`             | CADD SNV score file, for the `CADD` plugin. 81 GB — worth overriding with a local path. GRCh38 only                      |
| `--vep_cadd_snv_tbi`         | Index for `--vep_cadd_snv`. Required whenever `--vep_cadd_snv` is set                                                    |
| `--vep_cadd_indel`           | CADD indel score file, for the `CADD` plugin. GRCh38 only                                                                |
| `--vep_cadd_indel_tbi`       | Index for `--vep_cadd_indel`. Required whenever `--vep_cadd_indel` is set                                                |
| `--vep_revel`                | REVEL release zip, or a prepared score file, for the `REVEL` plugin. GRCh38 only                                         |
| `--vep_revel_tbi`            | Index for `--vep_revel`. Required only when `--vep_revel` is an already-prepared file                                    |
| `--vep_eve`                  | EVE release zip, or a merged VCF, for the `EVE` plugin. **Not** on by default. GRCh38 only. Default = `null`             |
| `--vep_eve_tbi`              | Index for `--vep_eve`. Required only when `--vep_eve` is an already-merged VCF. Default = `null`                         |

#### Minimap2 Options

| Parameter                    | Description                                                                                 |
| ---------------------------- | ------------------------------------------------------------------------------------------- |
| `--minimap2_ont_model`       | specifies which model to use minimap2 with for ONT samples. Default = `null`                |
| `--minimap2_pb_model`        | specifies which model to use minimap2 with for PacBio samples. Default = `null`             |
| `--save_secondary_alignment` | A boolean to specify if secondary alignments are kept in aligned bam file. Default = `true` |

#### ASCAT Options

| Parameter                     | Description                                                                                                                                                                                                                 |
| ----------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--ascat_ploidy`              | integer to enforce a given ploidy value. Default = `null`                                                                                                                                                                   |
| `--ascat_purity`              | integer to enforce a given purity value. Default = `null`                                                                                                                                                                   |
| `--ascat_min_base_qual`       | integer to specify a minimum base quality for ascat's allele counter. Default = `20`                                                                                                                                        |
| `--ascat_min_counts`          | integer to specify a minimum number of counts for ascat's allele counter. Default = `10`                                                                                                                                    |
| `--ascat_min_map_qual`        | integer to specify a minimum mapping quality for ascat's allele counter. Default = `10`                                                                                                                                     |
| `--ascat_penalty`             | integer to specify a penalty value for ascat. Default = `150`                                                                                                                                                               |
| `--ascat_longread_bins`       | integer to specify the binsize for ascat long reads. Default = `2000`                                                                                                                                                       |
| `--ascat_allelecounter_flags` | flags to pass to ascat's allele counter. Default = `"-f 0"`                                                                                                                                                                 |
| `--ascat_chroms`              | string to enforce a subset of chromosomes on the sample, ie `"(c(1:21,'X','Y')). Default = `null`                                                                                                                           |
| `--ascat_allele_files`        | A full path to a zipped folder containing allele files for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Must be zipped and have `.zip` format. Default = `null`                             |
| `--ascat_loci_files`          | A full path to a zipped folder containing loci files for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Must be zipped and have `.zip` format. Default = `null`                               |
| `--ascat_gc_file`             | A full path to a GC correction file for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Optionally can be zipped and have either `.txt` or `.txt.zip` format. Default = `null`                 |
| `--ascat_rt_file`             | A full path to a replication timing correction file for [ASCAT](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WGS). Optionally can be zipped and have either `.txt` or `.txt.zip` format. Default = `null` |
| `--ascat_pdf_plots`           | string to enable output pltos in pdf format. Default = `false`                                                                                                                                                              |

#### Fibertools Options

| Parameter           | Description                                                                      |
| ------------------- | -------------------------------------------------------------------------------- |
| `--autocorrelation` | A boolean to enable autocorrelation computation in fibertools. Default = `false` |

#### SEVERUS Options

| Parameter              | Description                                                                          |
| ---------------------- | ------------------------------------------------------------------------------------ |
| `--severus_minsupport` | Minimum number of supporting reads required for SEVERUS to call an SV. Default = `3` |

#### Report Options

| Parameter             | Description                                                                                                                                                                                                                                                                                                                                |
| --------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--report_src`        | Override the report tool source tree (bin/, R/, templates/, assets/). Not needed for normal runs: a copy of [lrsomatic_report](https://github.com/ljwharbers/lrsomatic_report) ships inside the pipeline. Point it at a local checkout to render with an unreleased version of the tool. Default = `${projectDir}/assets/lrsomatic_report` |
| `--report_gene_panel` | Gene panel(s) applied when the report opens, as a comma-separated list. Each entry is `none` (no filtering), a builtin panel name (`lymphoid` or `sarcoma`), or a path to a TSV file with a `gene` column. Default = `null`, i.e. unfiltered                                                                                               |

Gene panel filtering is a view, not a filter on the data: every builtin panel is embedded in
the rendered report and the reader can tick and untick them (or clear them all for the
unfiltered table) in the browser. `--report_gene_panel` only decides which ones are ticked on
load. A custom panel is a tab-separated file with a header row containing at least a `gene`
column:

```tsv
gene	panel	note
TP53	mypanel	Tumour suppressor
KRAS	mypanel	Oncogene
```

A panel may also carry `chrom`, `start` and `end` columns — all three or none. With
coordinates, structural variants are matched on position (within 1 Mb of a breakend, or
100 kb of the SV span) rather than on the VEP gene symbol, which is what makes breakend
filtering reliable: whether a breakend carries a gene symbol at all depends on the VEP
invocation. A coordinate-carrying panel must declare the reference its coordinates are
valid for, either as a leading `# reference: hg38` comment or as a `reference` column; a
panel declaring a reference other than the one the sample was called against is a hard
error rather than a silently wrong filter. Symbol-only panels need no declaration. The
builtin panels ship one file per reference and are selected by their bare name
(`lymphoid`, `sarcoma`), resolved against the detected reference.

```bash
nextflow run IntGenomicsLab/lrsomatic \
    -profile <docker/singularity> \
    --input samplesheet.csv \
    --outdir results \
    --report_gene_panel /path/to/mypanel.tsv
```

##### Applying several panels at once

Pass a comma-separated list to open the report with several panels applied. They are
**unioned**: a variant or SV is kept if it hits any of them. Builtin names and custom paths
can be mixed freely.

```bash
    --report_gene_panel lymphoid,sarcoma
    --report_gene_panel 'lymphoid,/path/to/mypanel.tsv'
```

With two or more panels active, each `panel_hit` entry in the SV table gains a trailing
`[panel]` naming which one matched — under a union that is all that distinguishes two hits
on the same gene. With a single panel the labels read exactly as they always have.

An entry is read as a panel **file** if it contains a `/` or ends in `.tsv`, and as a
builtin panel name otherwise. In practice that means a custom panel needs a path or a
`.tsv` name — `--report_gene_panel mypanel` is looked up as a builtin even if a file called
`mypanel` sits next to you.

Three things are checked before the run starts, so a mistake costs seconds rather than a
full pipeline:

- `none` means unfiltered and cannot be combined with a real panel.
- A panel file that does not exist, and a builtin name that is not one of the bundled
  panels, are both errors — a typo cannot quietly produce an unfiltered report.
- Two panel files sharing a base name cannot be combined, whatever directories they live
  in — they are staged side by side and would collide. Rename one.

#### WAKHAN Options

| Parameter         | Description                                                                                             |
| ----------------- | ------------------------------------------------------------------------------------------------------- |
| `--wakhan_chroms` | A string specifying a subset of chromosomes for WAKHAN to process, e.g. `"chr1,chr2"`. Default = `null` |

#### Mutational Signature Options

Mutational signature analysis runs [SigProfilerMatrixGenerator](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator) on the PASS SNVs and indels of the phased somatic VCF of every sample (SBS, DBS and ID matrices at all context sizes, with plots) and then fits COSMIC reference signatures per sample with [SigProfilerAssignment](https://github.com/SigProfilerSuite/SigProfilerAssignment) (SBS96, DBS78 and ID83). It needs SigProfilerMatrixGenerator's per-genome payload (~3 GB, the transcriptional-strand-annotated chromosomes), which is not shipped with the pipeline. Either:

- run once with `--download_sigprofiler_genome`; the payload is installed, checksum-verified and published to `<outdir>/cache/sigprofiler/volume`, or
- point `--sigprofiler_genome_dir` at an existing SigProfilerMatrixGenerator volume (a directory containing `tsb/<genome>/`, e.g. one created with `SigProfilerMatrixGenerator install GRCh38 --volume <dir>` or the published cache from a previous run).

Running with neither (and without `--skip_signatures`) stops the pipeline at start-up with an explanatory error.

| Parameter                                   | Description                                                                                                                                                                                         |
| ------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--sigprofiler_genome_dir`                  | Full path to a SigProfilerMatrixGenerator volume containing `tsb/<sigprofiler_genome>/`. Default = `null`                                                                                           |
| `--download_sigprofiler_genome`             | A boolean to download and install the genome payload during the run (published to `<outdir>/cache/sigprofiler/volume`). Default = `false`                                                           |
| `--sigprofiler_cosmic_version`              | COSMIC reference signature version fitted by SigProfilerAssignment. Default = `3.6`                                                                                                                 |
| `--sigprofiler_exclude_signature_subgroups` | Comma-separated SigProfilerAssignment signature subgroups to exclude from the fit, e.g. `"Artifact_signatures,Lymphoid_signatures"` (see the SigProfilerAssignment documentation). Default = `null` |
| `--sigprofiler_matrix_args`                 | Extra arguments for `SigProfilerMatrixGenerator matrix_generator`. Default = `"--plot"`                                                                                                             |
| `--sigprofiler_assignment_args`             | Extra arguments for `SigProfilerAssignment cosmic_fit`, e.g. `"--make_plots False"`. Default = `null`                                                                                               |

The tools run from a purpose-built image (`ghcr.io/ljwharbers/sigprofiler`) because CHM13 support is not in a SigProfiler release yet: SigProfilerMatrixGenerator comes from the branch behind [SigProfilerSuite/SigProfilerMatrixGenerator#250](https://github.com/SigProfilerSuite/SigProfilerMatrixGenerator/pull/250) (adds the `CHM13-T2T` genome) and SigProfilerAssignment from [ljwharbers/SigProfilerAssignment](https://github.com/ljwharbers/SigProfilerAssignment/tree/chm13-t2t-support), which adds COSMIC SBS/DBS signatures renormalised to the CHM13 trinucleotide/dinucleotide composition (stock SigProfilerAssignment silently falls back to GRCh37 signatures for CHM13). Indel (ID83) signatures are not genome-normalised by COSMIC and always use the GRCh37 set. Conda is not supported for this step.

#### Variant Filtering and Combining Options

These options control how variants from multiple callers are filtered and merged.

| Parameter                      | Description                                                                                         |
| ------------------------------ | --------------------------------------------------------------------------------------------------- |
| `--germline_var_keep`          | Expression or threshold for retaining germline variants after calling. Default = `null`             |
| `--somatic_var_keep`           | Expression or threshold for retaining somatic variants after calling. Default = `null`              |
| `--germline_var_combine`       | Strategy for combining germline variant caller outputs (e.g. union, intersection). Default = `null` |
| `--somatic_var_combine`        | Strategy for combining somatic variant caller outputs (e.g. union, intersection). Default = `null`  |
| `--prioritize_caller_germline` | Comma-separated caller priority order used when combining germline calls. Default = `null`          |
| `--prioritize_caller_somatic`  | Comma-separated caller priority order used when combining somatic calls. Default = `null`           |

#### PON Options

| Parameter                | Description                                                                                                                                                                                                                                  |
| ------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--clairsto_pon_vcfs`    | Full path to one or more Panel of Normals VCF files for ClairS-TO small variant filtering. Default = `null`                                                                                                                                  |
| `--clairsto_pon_flags`   | Population allele matching flags for ClairS-TO PON VCFs (one per VCF, comma-separated). Default = `null`                                                                                                                                     |
| `--deepsomatic_pon_vcfs` | Full path to one or more bgzipped, tabix-indexed PON VCF files (for example, `.vcf.gz`) passed to DeepSomatic `--population_vcfs`. If not set, uses container-bundled defaults in tumor-only mode or no PON in paired mode. Default = `null` |

#### Advanced Options

| Parameter         | Description                                                                                                                                                                 |
| ----------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--use_gpu`       | A boolean to enable GPU acceleration for DeepVariant and DeepSomatic. Requires a GPU-enabled compute environment. Default = `false`                                         |
| `--generate_gvcf` | A boolean to enable gVCF output from DeepVariant (germline) and DeepSomatic (somatic). gVCF files include calls at all positions, not just variant sites. Default = `false` |

#### Genome-Derived Parameters

The following parameters are automatically populated from the `--genome` iGenomes configuration and do not normally need to be set manually. They can be overridden when using a custom genome or reference build not present in the iGenomes configuration.

| Parameter                  | Description                                                                                                                                                                   |
| -------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--fasta`                  | Full path to the reference FASTA file. Auto-populated from `--genome`. Override for custom genomes.                                                                           |
| `--bed_file`               | BED file of callable/target regions passed to SEVERUS for SV calling. Auto-populated from `--genome`.                                                                         |
| `--pon_file`               | Panel of Normals breakpoint table (bgzipped CSV) for SEVERUS somatic SV filtering in tumor-only mode. Auto-populated from `--genome`.                                         |
| `--centromere_bed`         | BED file of centromere coordinates passed to WAKHAN. Auto-populated from `--genome`.                                                                                          |
| `--genome_name`            | Assembly name string passed to ASCAT for genome-specific reference file selection. Auto-populated from `--genome`.                                                            |
| `--vep_genome`             | VEP genome identifier (e.g. `GRCh38`, `T2T-CHM13v2.0`). Auto-populated from `--genome`. Override for CHM13 or custom assemblies.                                              |
| `--vep_species`            | VEP species identifier. Auto-populated from `--genome`. Override for non-standard assemblies (e.g. `homo_sapiens_gca009914755v4` for CHM13).                                  |
| `--sigprofiler_genome`     | SigProfilerMatrixGenerator / SigProfilerAssignment genome name (`GRCh38` or `CHM13-T2T`). Auto-populated from `--genome`.                                                     |
| `--sigprofiler_genome_url` | URL of the `<sigprofiler_genome>.tar.gz` payload for `--download_sigprofiler_genome`; unset means the AlexandrovLab FTP. Auto-populated from `--genome` (set for CHM13 only). |

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull IntGenomicsLab/lrsomatic
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [IntGenomicsLab/lrsomatic releases page](https://github.com/IntGenomicsLab/lrsomatic/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## VEP plugins

VEP is given extra pathogenicity and clinical-significance annotation through plugins. On
`--genome GRCh38` and `--genome CHM13` these are **on by default**: the pipeline resolves the right
resources for the assembly, downloads them, and reshapes the ones that cannot be handed to VEP as
published. Nothing has to be prepared by hand first.

To annotate with VEP alone, pass `--skip_vep_plugins`. The plugins are also off whenever no default
resolves — under `--igenomes_ignore`, or with a `--genome` other than GRCh38 and CHM13.

The pipeline hosts and redistributes none of this data: every default points at the resource's
original source. Complying with each licence remains your responsibility, and **CADD, REVEL and EVE
are free for non-commercial use only** — see the callout at the end of this section.

Plugins are applied to the germline and somatic VEP runs. They are deliberately **not** applied to
the structural-variant run, since missense and splice scores carry no meaning on SEVERUS breakends.

### What is enabled by default, per assembly

| Tool              | GRCh38                                   | CHM13                         | Default download                       |
| ----------------- | ---------------------------------------- | ----------------------------- | -------------------------------------- |
| **SIFT**          | already in the VEP cache, no file needed | `PolyPhen_SIFT` plugin        | CHM13: 13 GB database                  |
| **PolyPhen**      | already in the VEP cache, no file needed | `PolyPhen_SIFT` plugin        | as above, the same database            |
| **AlphaMissense** | `AlphaMissense` plugin                   | `AlphaMissenseProtein` plugin | 613 MB (GRCh38) / 1.1 GB (CHM13)       |
| **ClinVar**       | `--custom` annotation                    | `--custom`, CHM13-lifted VCF  | 105 MB (GRCh38) / 190 MB (CHM13)       |
| **CADD**          | `CADD` plugin                            | not available — see below     | 81 GB SNVs + 1.2 GB indels             |
| **REVEL**         | `REVEL` plugin                           | not available — see below     | 667 MB release zip                     |
| **EVE**           | `EVE` plugin, opt-in                     | not available — see below     | none — `--vep_eve` enables it (9.6 GB) |

That comes to roughly **83.5 GB on GRCh38** and **14.3 GB on CHM13** for a first run, almost all of
it CADD and the pangenome database respectively. On GRCh38, SIFT and PolyPhen already come out of
the VEP cache through the default `--everything`, so `--vep_polyphen_sift_db` is only used on CHM13.

EVE is the one resource left opt-in: its release is a 9.6 GB zip of per-protein VCFs that then has
to be merged, which is a lot of work for a predictor AlphaMissense largely covers. Enable it with
`--vep_eve https://evemodel.org/api/proteins/bulk/download/` (or a path to the zip, or to a VCF you
have already merged).

### Downloading and reusing the data

Prepared files are **not** published to `--outdir`. They live in the Nextflow work directory, where
they are shared by every VEP task in the run and reused by a `-resume` against the same work
directory — but a fresh work directory downloads and prepares them again. Only REVEL and EVE are
prepared at all; the two AlphaMissense tables are fetched already indexed, which spares the
reshaping rather than the download, since the table is about the size of the release it replaces.

For repeat runs, download once and pass local paths, which skips both the download and the
preparation step:

```bash
# once
bin/prepare_vep_plugin_data.sh revel /data/vep/revel-v1.3_all_chromosomes.zip /data/vep

# then, on every run
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv --outdir ./results --genome GRCh38 \
  --vep_cadd_snv /data/vep/whole_genome_SNVs.tsv.gz \
  --vep_cadd_snv_tbi /data/vep/whole_genome_SNVs.tsv.gz.tbi \
  --vep_revel /data/vep/revel_grch38.tsv.gz \
  --vep_revel_tbi /data/vep/revel_grch38.tsv.gz.tbi \
  -profile singularity
```

Everything except REVEL and EVE needs no preparation — download it with `curl` and pass the path.
`bin/prepare_vep_plugin_data.sh` is what the pipeline itself runs for those two, so a file you build
with it is identical to one it would prepare; run it with `--help` for usage.

Whether an index is required depends on the shape of what you supply:

- **AlphaMissense, ClinVar and CADD** are used exactly as given, so their `_tbi` parameter is
  always required alongside them. Overriding a data file **drops the default index**, deliberately:
  an index belongs to the file it was built from, and pairing yours with ours would point tabix at
  an index for different content. Supply both, or neither.
- **REVEL and EVE** ship as zip archives. Pass a `.zip` and the pipeline unpacks and reshapes it;
  pass a prepared file and its index instead to use it directly. A remote `.zip` is fetched with
  `wget` rather than staged like the other files, because neither release host can be staged by
  Nextflow directly: REVEL's answers `403` to a request carrying no `User-Agent`, and EVE's
  redirects HTTPS to HTTP, which Nextflow refuses to follow. EVE's 9.6 GB archive also serves at a
  few hundred KB/s, so expect hours — another reason it is opt-in.

### Where each default comes from

| Resource                | Source                                                                                                                            | Size   | Prepared by the pipeline                   | Licence                      |
| ----------------------- | --------------------------------------------------------------------------------------------------------------------------------- | ------ | ------------------------------------------ | ---------------------------- |
| AlphaMissense (GRCh38)  | `https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz`, with an index we host                                | 613 MB | no, the index is fetched ready             | CC BY 4.0                    |
| AlphaMissense (protein) | a gene-symbol-keyed table we host, built from the AlphaMissense protein-space release                                             | 1.1 GB | no, fetched ready — see `CITATIONS.md`     | CC BY 4.0                    |
| Pangenome PolyPhen/SIFT | `https://ftp.ensembl.org/pub/current_variation/pangenomes/Human/homo_sapiens_pangenome_PolyPhen_SIFT_20240502.db`                 | 13 GB  | no, it is an SQLite database               | Ensembl / EMBL-EBI open      |
| ClinVar (GRCh38)        | `https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz` (+ `.tbi`)                                                   | 105 MB | no, `.tbi` is published                    | public domain                |
| ClinVar (CHM13)         | `clinvar_20240624_GCA_009914755.4.vcf.gz` under `https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/` | 190 MB | no, `.tbi` is published                    | public domain                |
| CADD v1.7 SNVs          | `https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz` (+ `.tbi`)                                 | 81 GB  | no, `.tbi` is published                    | free for non-commercial use  |
| CADD v1.7 indels        | `https://krishna.gs.washington.edu/download/CADD/v1.7/GRCh38/gnomad.genomes.r4.0.indel.tsv.gz` (+ `.tbi`)                         | 1.2 GB | no, `.tbi` is published                    | free for non-commercial use  |
| REVEL v1.3              | `https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip`                                                                | 667 MB | unpacked, re-sorted on GRCh38, indexed     | free for non-commercial use  |
| EVE (opt-in)            | `https://evemodel.org/api/proteins/bulk/download/`                                                                                | 9.6 GB | unpacked, per-protein VCFs merged, indexed | non-commercial; see the site |

### Example: GRCh38

Defaults do the work, so nothing plugin-related is needed:

```bash
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv \
  --outdir ./results \
  --genome GRCh38 \
  -profile docker
```

### Example: CHM13

```bash
nextflow run IntGenomicsLab/lrsomatic \
  --input samplesheet.csv \
  --outdir ./results \
  --genome CHM13 \
  --vep_cache $HOME/.vep \
  --vep_cache_version 107 \
  -profile docker
```

Passing a GRCh38-only resource together with a CHM13 cache is rejected before the run starts,
rather than failing hours later inside VEP or after a long download.

### How annotation reaches CHM13 at all

The CHM13 cache is the Ensembl rapid-release HPRC cache: gene and transcript models only, with no
variation, RefSeq or regulatory data. Every genome-coordinate-keyed score resource is published for
GRCh37/GRCh38 alone, and there are no lifted score tracks for the assembly elsewhere either.

Two of these predictors get there anyway, because they score _proteins_ rather than genome positions:

- **`PolyPhen_SIFT`** keys its lookup on the MD5 of the peptide sequence, which is exactly why
  Ensembl ships a pangenome database covering the HPRC assemblies.
- **`AlphaMissenseProtein`** (in `assets/vep_plugins/`) keys on gene symbol plus amino-acid
  substitution, using a table built from AlphaMissense's own protein-space release. A row is used
  only when both the reference and the alternate amino acid match what VEP computed for the CHM13
  transcript; where the proteins genuinely disagree it reports `aa_mismatch` and returns no score
  rather than a score for the wrong substitution. Check the `AlphaMissenseProtein_match` field to
  see how each lookup resolved.

### What is not available, and why

- **MutationTaster** — there is no MutationTaster plugin in Ensembl's `VEP_plugins`, and
  MutationTaster 2021 is a web service. Its scores are redistributed through dbNSFP, so an offline
  route does exist, but dbNSFP's academic-only terms and the separately licensed components it
  bundles make it unsuitable as a pipeline default.
- **CADD on CHM13** — CADD scores non-coding positions as well as coding ones, so unlike SIFT,
  PolyPhen and AlphaMissense it has no protein-space representation to fall back on. It is
  structurally unavailable on CHM13, not merely unpublished.
- **REVEL and EVE on CHM13** — both could in principle be re-keyed into protein space, but neither
  publishes licence terms that clearly permit distributing a derived table. AlphaMissense covers the
  same class of variant and is CC BY 4.0, so it is used instead.
- **SpliceAI** — not currently wired up on either assembly.

> [!IMPORTANT]
> **CADD and REVEL are enabled by default and are free for non-commercial use only**; EVE, if you
> enable it, is the same. The pipeline cannot accept those terms on your behalf: if your work is
> commercial, pass `--skip_vep_plugins`, or set only the resources you are licensed for.

> [!IMPORTANT]
> Check that the contig naming of every file you supply matches your reference. The GRCh38 reference
> used here is GATK-style (`chr1`), while NCBI's ClinVar VCF ships Ensembl-style names (`1`). A
> mismatch produces empty annotation rather than an error, so it is easy to miss — confirm that
> `CLNSIG` values actually appear in an annotated VCF before trusting them.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow `24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
