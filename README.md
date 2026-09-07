# IntGenomicsLab/lrsomatic

[![Open in GitHub Codespaces](https://img.shields.io/badge/Open_In_GitHub_Codespaces-black?labelColor=grey&logo=github)](https://github.com/codespaces/new/IntGenomicsLab/lrsomatic)
[![GitHub Actions CI Status](https://github.com/IntGenomicsLab/lrsomatic/actions/workflows/nf-test.yml/badge.svg)](https://github.com/IntGenomicsLab/lrsomatic/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/IntGenomicsLab/lrsomatic/actions/workflows/linting.yml/badge.svg)](https://github.com/IntGenomicsLab/lrsomatic/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.17751829-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.17751829)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.04.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/IntGenomicsLab/lrsomatic)

## Introduction

**IntGenomicsLab/lrsomatic** is a robust bioinformatics pipeline designed for processing and analyzing **somatic DNA sequencing** data for long-read sequencing technologies from **Oxford Nanopore** and **PacBio**. It supports both canonical base DNA and modified base calling, including specialized applications such as **Fiber-seq**.

This **end-to-end pipeline** handles the entire workflow — **from raw read processing and alignment, to comprehensive somatic variant calling**, including single nucleotide variants, indels, structural variants, copy number alterations, and modified bases.

It can be run in both **matched tumour-normal** and **tumour-only mode**, offering flexibility depending on the users study design.

Developed using **Nextflow DSL2**, it offers high portability and scalability across diverse computing environments. By leveraging Docker or Singularity containers, installation is streamlined and results are highly reproducible. Each process runs in an isolated container, simplifying dependency management and updates. Where applicable, pipeline components are sourced from **nf-core/modules**, promoting reuse, interoperability, and consistency within the broader Nextflow and nf-core ecosystems.

For more information on how to run the pipeline, you can also go [here](https://intgenomicslab.github.io/lrsomatic).

## Pipeline summary

![image](./assets/lrsomatic_1.0.png)

**1) Pre-processing:**

a. Raw read QC ([`cramino`](https://github.com/wdecoster/cramino))

b. Alignment to the reference genome ([`minimap2`](https://github.com/lh3/minimap2))

c. Post alignment QC ([`cramino`](https://github.com/wdecoster/cramino), [`samtools idxstats`](https://github.com/samtools/samtools), [`samtools flagstats`](https://github.com/samtools/samtools), [`samtools stats`](https://github.com/samtools/samtools))

d. Specific for calling modified base calling ([`Modkit`](https://github.com/nanoporetech/modkit), [`Fibertools`](https://github.com/fiberseq/fibertools-rs))

**2i) Matched mode: small variant calling:**

a. Calling Germline SNPs ([`Clair3`](https://github.com/HKU-BAL/Clair3))

b. Phasing and Haplotagging the SNPs in the normal and tumour BAM ([`LongPhase`](https://github.com/twolinin/longphase))

c. Calling somatic SNVs ([`ClairS`](https://github.com/HKU-BAL/ClairS))

**2ii) Tumour only mode: small variant calling:**

a. Calling Germline SNPs and somatic SNVs ([`ClairS-TO`](https://github.com/HKU-BAL/ClairS-TO))

b. Phasing and Haplotagging germline SNPs in tumour BAM ([`LongPhase`](https://github.com/twolinin/longphase))

**3) Large variant calling:**

a. Somatic structural variant calling ([`Severus`](https://github.com/KolmogorovLab/Severus))

b. Somatic structural variant and copy-number calling ([`SAVANA`](https://github.com/cortes-ciriano-lab/savana))

c. Copy number alterion calling; long read version of ([`ASCAT`](https://github.com/VanLoo-lab/ascat))

**4) Annotation:**

a. Small variant annotation ([`VEP`](https://github.com/Ensembl/ensembl-vep))

b. Structural variant annotation ([`VEP`](https://github.com/Ensembl/ensembl-vep))

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/guidelines/graphic_design/workflow_diagrams#examples for examples.   -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First prepare a samplesheet with your input data that looks as follows:

```csv
sample,bam_tumor,bam_normal,platform,sex,fiber
sample1,tumour.bam,normal.bam,ont,female,n
sample2,tumour.bam,,ont,female,y
sample3,tumour.bam,,pb,male,n
sample4,tumour.bam,normal.bam,pb,male,y
```

Each row represents a sample. The bam files should always be unaligned bam files. All fields except for `bam_normal` are required. If `bam_normal` is empty, the pipeline will run in tumour only mode. `platform` should be either `ont` or `pb` for Oxford Nanopore Sequencing or PacBio sequencing, respectively. `sex` refers to the biological sex of the sample and should be either `female` or `male`. Finally, `fiber` specifies whether your sample is Fiber-seq data or not and should have either `y` for Yes or `n` for No.

Now, you can run the pipeline using:

```bash
nextflow run IntGenomicsLab/lrsomatic \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

More detail is given in our [usage documentation](/docs/usage.md)

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

IntGenomicsLab/lr_somatic was originally written by Luuk Harbers, Robert Forsyth, Alexandra Pančíková, Marios Eftychiou, Ruben Cools, Laurens Lambrechts, and Jonas Demeulemeester.

## Pipeline output

This pipeline produces a series of different output files. The main output is an aligned and phased tumour bam file. This bam file can be used by any typical downstream tool that uses bam files as input. Furthermore, we have sample-specific QC outputs from `cramino` (fastq), `cramino` (bam), `mosdepth`, `samtools` (stats/flagstat/idxstats), and optionally `fibertools`. Finally, we have a `multiqc` report from that combines the output from `mosdepth` and `samtools` into one html report.

Besides QC and the aligned and phased bam file, we have output from (structural) variant and copy number callers, of which some are optional. The output from these variant callers can be found in their respective folders. For small and structural variant callers (`clairS`, `clairS-TO`, `severus`, and `savana`) these will contain, among others, `vcf` files with called variants. For `ascat` and `savana` these also contain files with copy number information.

Example output directory structure:

```
├── Sample 1
│    ├── ascat
│    ├── bamfiles
│    ├── qc
│    │    ├── tumor
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── samtools
│    ├── variants
│    │   ├──clairS-TO
│    │   ├──severus
│    │   ├──savana
│    ├── vep
│    │   ├── germline
│    │   ├── somatic
│    │   ├── SVs
│
├── Sample 2
│    ├── ascat
│    ├── bamfiles
│    ├── qc
│    │    ├── tumor
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── samtools
│    │    ├── normal
│    │    │   ├── cramino_aln
│    │    │   ├── cramino_ubam
│    │    │   ├── fibertoolsrs
│    │    │   ├── mosdepth
│    │    │   ├── samtools
│    ├── variants
│    │   ├── clair3
│    │   ├── clairS
│    │   ├── severus
│    │   ├── savana
│    ├── vep
│    │   ├── germline
│    │   ├── somatic
│    │   ├── SVs
├── pipeline_info
```

more detail is given in our [output documentation](/docs/output.md)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

If you use `IntGenomicsLab/lrsomatic` for your analysis, please cite it using the following:

> LRSomatic: a highly scalable and robust pipeline for somatic variant calling in long-read sequencing data
>
> Robert A. Forsyth*, Luuk Harbers*, Amber Verhasselt, Ana-Lucía Rocha Iraizós, Sidi Yang, Joris Vande Velde, Christopher Davies, Nischalan Pillay, Laurens Lambrechts, Jonas Demeulemeester
>
> bioRxiv 2026.02.26.707772; doi: https://doi.org/10.64898/2026.02.26.707772

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
