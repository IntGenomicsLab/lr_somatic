/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_lrsomatic_pipeline'
include { getGenomeAttribute     } from '../subworkflows/local/utils_nfcore_lrsomatic_pipeline'

//
// IMPORT MODULES
//
include { SAMTOOLS_MERGE                    } from '../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGE } from '../modules/nf-core/samtools/index/main'
include { MINIMAP2_INDEX                    } from '../modules/nf-core/minimap2/index/main'
include { MINIMAP2_ALIGN                    } from '../modules/nf-core/minimap2/align/main'
include { CRAMINO as CRAMINO_PRE            } from '../modules/local/cramino/main'
include { CRAMINO as CRAMINO_POST           } from '../modules/local/cramino/main'
include { NANOPLOT as NANOPLOT_PRE          } from '../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_POST         } from '../modules/nf-core/nanoplot/main'
include { MOSDEPTH                          } from '../modules/nf-core/mosdepth/main'
include { ASCAT                             } from '../modules/nf-core/ascat/main'
include { SEVERUS                           } from '../modules/nf-core/severus/main.nf'
include { METAEXTRACT                       } from '../modules/local/metaextract/main'
include { WAKHAN                            } from '../modules/local/wakhan/main'
include { FIBERTOOLSRS_PREDICTM6A           } from '../modules/local/fibertoolsrs/predictm6a'
include { FIBERTOOLSRS_FIRE                 } from '../modules/local/fibertoolsrs/fire'
include { FIBERTOOLSRS_NUCLEOSOMES          } from '../modules/local/fibertoolsrs/nucleosomes'
include { FIBERTOOLSRS_QC                   } from '../modules/local/fibertoolsrs/qc'
include { ENSEMBLVEP_VEP as SOMATIC_VEP     } from '../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as GERMLINE_VEP    } from '../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as SV_VEP          } from '../modules/nf-core/ensemblvep/vep/main.nf'
include { ENSEMBLVEP_VEP as VEP_SAVANA      } from '../modules/nf-core/ensemblvep/vep/main.nf'
include { WHATSHAP_STATS                    } from '../modules/nf-core/whatshap/stats/main'
include { MODKIT_PILEUP                     } from '../modules/nf-core/modkit/pileup/main'

//
// IMPORT SUBWORKFLOWS
//
include { PREPARE_REFERENCE_FILES         } from '../subworkflows/local/prepare_reference_files'
include { PREPARE_ANNOTATION              } from '../subworkflows/local/prepare_annotation'
include { BAM_STATS_SAMTOOLS              } from '../subworkflows/nf-core/bam_stats_samtools/main'
include { TUMORONLY_SMALLVAR              } from '../subworkflows/local/tumor_only/tumoronly_smallvar'
include { PAIRED_SMALLVAR_SOMATIC         } from '../subworkflows/local/paired/paired_smallvar_somatic'
include { PAIRED_SMALLVAR_GERMLINE        } from '../subworkflows/local/paired/paired_smallvar_germline'
include { PHASING_HAPLOTYPING             } from '../subworkflows/local/phasing_haplotyping'
include { TUMORONLY_SAVANA                } from '../subworkflows/local/tumor_only/tumoronly_savana'
include { PAIRED_SAVANA                   } from '../subworkflows/local/paired/paired_savana'



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LRSOMATIC {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    // Channel format is [[meta], [bam]].
    // Where [meta] is [id, paired_data, method, specs, type]

    main:

    def clair3_modelMap = [
        'dna_r10.4.1_e8.2_400bps_sup@v5.2.0': 'r1041_e82_400bps_sup_v520',
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'r1041_e82_400bps_sup_v500',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'r1041_e82_400bps_sup_v430',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'r1041_e82_400bps_sup_v420',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'r1041_e82_400bps_sup_v410',
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'r1041_e82_260bps_sup_v400',
        'hifi_revio'                        : 'hifi'
    ]

    def clairs_modelMap = [
        'dna_r10.4.1_e8.2_260bps_sup@v4.0.0': 'ont_r10_dorado_sup_4khz',
        'dna_r10.4.1_e8.2_400bps_sup@v4.1.0': 'ont_r10_dorado_sup_4khz',
        'dna_r10.4.1_e8.2_400bps_sup@v4.2.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v4.3.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v5.0.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'dna_r10.4.1_e8.2_400bps_sup@v5.2.0': 'ont_r10_dorado_sup_5khz_ssrs',
        'hifi_revio'                        : 'hifi_revio_ssrs'

    ]

    // Load in igenomes
    params.fasta = getGenomeAttribute('fasta')
    params.genome_name = getGenomeAttribute('genome_name')
    params.ascat_allele_files = getGenomeAttribute('ascat_alleles')
    params.ascat_loci_files = getGenomeAttribute('ascat_loci')
    params.ascat_gc_file = getGenomeAttribute('ascat_loci_gc')
    params.ascat_rt_file = getGenomeAttribute('ascat_loci_rt')
    params.centromere_bed = getGenomeAttribute('centromere_bed')
    params.pon_file = getGenomeAttribute('pon_file')
    params.bed_file = getGenomeAttribute('bed_file')
    params.savana_contigs = getGenomeAttribute('savana_contigs')
    params.savana_g1000_vcf = getGenomeAttribute('savana_g1000_vcf')
    params.vep_genome = getGenomeAttribute('vep_genome')
    params.vep_species = getGenomeAttribute('vep_species')

    // Convert comma-separated caller strings to lists for internal use
    params.germline_var_keep = params.germline_var_keep instanceof List
        ? params.germline_var_keep
        : params.germline_var_keep.tokenize(',').collect { it.trim() }
    params.somatic_var_keep = params.somatic_var_keep instanceof List
        ? params.somatic_var_keep
        : params.somatic_var_keep.tokenize(',').collect { it.trim() }

    if (params.clairsto_pon_vcfs != null) {
        pon_files = params.clairsto_pon_vcfs.split(',').collect { f -> file(f.trim()) }
        if (params.clairsto_pon_flags != null) {
            pon_flags = params.clairsto_pon_flags.split(',').collect { f -> f.trim() }
        } else if (params.genome == 'GRCh38') {
            pon_flags = ["True", "True", "False", "False"]
        } else if (params.genome == 'CHM13') {
            pon_flags = ["True", "True", "False", "False", "False"]
        } else {
            pon_flags = pon_files.collect { "False" }
        }
    }
    else if (params.genome == 'GRCh38') {
        pon_files  = [
            getGenomeAttribute('gnomad'),
            getGenomeAttribute('dbsnp'),
            getGenomeAttribute('onekgenomes'),
            getGenomeAttribute('colors'),
        ]
        pon_flags = [
            "True",
            "True",
            "False",
            "False"
        ]
    }
    else if (params.genome == 'CHM13') {
        pon_files  = [
            getGenomeAttribute('gnomad'),
            getGenomeAttribute('dbsnp'),
            getGenomeAttribute('onekgenomes'),
            getGenomeAttribute('colors'),
            getGenomeAttribute('asap')
        ]
        pon_flags = [
            "True",
            "True",
            "False",
            "False",
            "False"
        ]
    }
    if (pon_files.size() != pon_flags.size()) {
        error "PoN VCFs and allele flags must have same length"
    }
    channel
        .of( tuple(pon_files, pon_flags) )
        .set { clairsto_pon_channel }
    // clairsto_pon_channel: [ [pon_vcf_path, ...], [is_population_allele_flag, ...] ]
    //   -- single tuple of parallel lists; each flag indicates whether the corresponding VCF
    //      is a population allele database (True) vs. a panel-of-normals artefact file (False)

    // DeepSomatic PON channel: user-supplied VCF paths, or empty list (process falls back to container defaults)
    ds_pon_files = params.deepsomatic_pon_vcfs != null
        ? params.deepsomatic_pon_vcfs.split(',').collect { f -> file(f.trim()) }
        : params.genome == 'CHM13'
            ? [
                getGenomeAttribute('gnomad'),
                getGenomeAttribute('dbsnp'),
                getGenomeAttribute('onekgenomes'),
                getGenomeAttribute('colors'),
                getGenomeAttribute('asap')
              ]
            : []
    // DeepSomatic requires no chromosome overlap across population VCFs.
    // When multiple databases are provided (e.g., CHM13 gnomad + 1kgenomes + colors + dbsnp + asap),
    // the merge is done inline inside DEEPSOMATIC_MAKEEXAMPLES and DEEPSOMATIC_POSTPROCESSVARIANTS
    // so that both callers can start in parallel as soon as BAMs are ready.
    channel.value( [[:], ds_pon_files] ).set { ds_pon_channel }
    // ds_pon_channel: [[:], [vcf_path, ...]] or [[:], []]
    //   -- raw unmerged PON VCF paths (no .tbi required); merging happens inline in each DeepSomatic process
    //   -- GRCh38/other + no user PON: empty list => process uses container-bundled GRCh38 defaults (tumor-only)

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

    //
    // MODULE: METAEXTRACT
    //
    // extracts the base calling model from the bam files

    // MODULE: METAEXTRACT (label: process_single)
    // Input:  [meta, [bam...]]
    METAEXTRACT( ch_samplesheet )

    basecall_meta = METAEXTRACT.out.meta_ext
    // basecall_meta: [meta, basecall_model_str, kinetics_str]
    //   basecall_model_str -- e.g. "dna_r10.4.1_e8.2_400bps_sup@v5.0.0" or "hifi_revio"
    //   kinetics_str       -- "true" if PacBio kinetics tags present, else "false"

    ch_samplesheet
        .join(basecall_meta)
        .map { meta, bam, basecall_model_meta, kinetics_meta ->
            def chosen_clair3_model = meta.clair3_model ?: clair3_modelMap.get(basecall_model_meta)
            def chosen_clairSTO_model = meta.clairSTO_model ?: clairs_modelMap.get(basecall_model_meta)
            def chosen_clairS_model = meta.clairS_model ?: clairs_modelMap.get(basecall_model_meta)
            def meta_new =[ id: meta.id,
                            paired_data: meta.paired_data,
                            type: meta.type,
                            platform: meta.platform,
                            sex: meta.sex,
                            fiber: meta.fiber,
                            replicate: meta.replicate,
                            n_replicates: meta.n_replicates,
                            clair3_model: chosen_clair3_model,
                            clairS_model: chosen_clairS_model,
                            clairSTO_model: chosen_clairSTO_model,
                            kinetics: kinetics_meta]
            return[ meta_new, bam ]
        }
        .set{ch_samplesheet}
    // ch_samplesheet (updated): [meta, [bam...]]
    //   meta fields: id, paired_data, type, platform, sex, fiber, replicate,
    //                clair3_model, clairS_model, clairSTO_model, kinetics
    //   bams are grouped per sample (multiple runs merged into a list)

    //
    // SUBWORKFLOW: PREPARE_REFERENCE_FILES
    // Decompresses the reference FASTA if needed, indexes it, downloads Clair3 models,
    // and decompresses ASCAT reference files
    // Input:  params.fasta, ASCAT file paths, basecall_meta, clair3_modelMap
    // Output: .prepped_fasta           -- [[:], fasta]
    //         .prepped_fai             -- [[:], fai]
    //         .downloaded_clair3_models-- [meta(id=model_name), model_dir]
    //         .allele_files / .loci_files / .gc_file / .rt_file  -- flat file collections
    //

    PREPARE_REFERENCE_FILES (
        params.fasta,
        params.ascat_allele_files,
        params.ascat_loci_files,
        params.ascat_gc_file,
        params.ascat_rt_file,
        basecall_meta,
        clair3_modelMap
    )

    downloaded_clair3_models = PREPARE_REFERENCE_FILES.out.downloaded_clair3_models
    // downloaded_clair3_models: [meta(id=clair3_model_name), model_dir]

    ch_nanoplot_pre_txt = channel.empty()

    if (!params.skip_qc && !params.skip_cramino) {

        //
        // MODULE: CRAMINO_PRE (label: process_medium)
        // Input:  [meta, [bam...]]  -- pre-alignment unaligned BAMs
        // Output: cramino_pre.out.arrow -- [meta, arrow_file] (feather format stats)
        //

        CRAMINO_PRE( ch_samplesheet )

        if (!params.skip_nanoplot) {

            //
            // MODULE: NANOPLOT_PRE (label: process_medium)
            // Input:  CRAMINO_PRE.out.arrow -- [meta, arrow_file]
            // Output: nanoplot HTML/txt reports
            //

            NANOPLOT_PRE(CRAMINO_PRE.out.arrow)

        }

    }

    // Each replicate is aligned separately so that minimap2 can tag its reads with a
    // unique @RG (sample + type + replicate).  Replicates are merged after alignment.
    // ch_samplesheet is [meta_with_replicate, [bam]] -- one item per replicate per sample.
    // ch_ubams defaults to ch_samplesheet; the fiber-seq block below may override it.
    ch_ubams = ch_samplesheet

    vep_cache = channel.empty()

    if (!params.skip_vep) {

        // SUBWORKFLOW: PREPARE_ANNOTATION
        // Validates or downloads the VEP cache directory
        // Output: .vep_cache -- path to VEP cache root directory
        PREPARE_ANNOTATION (
            params.vep_cache,
            params.vep_cache_version,
            params.vep_genome,
            params.vep_args,
            params.vep_species,
            params.download_vep_cache
        )
        ch_versions = ch_versions.mix(PREPARE_ANNOTATION.out.versions)
        // Wrap VEP cache path in a tuple with empty meta for use in ENSEMBLVEP_VEP
        vep_cache = PREPARE_ANNOTATION.out.vep_cache.map {cache -> [[:], cache] }
        // vep_cache: [[:], cache_dir_path]  -- empty meta + VEP cache directory

    }

    ch_versions = ch_versions.mix(PREPARE_REFERENCE_FILES.out.versions)
    ch_fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta  // [[:], fasta]
    ch_fai   = PREPARE_REFERENCE_FILES.out.prepped_fai    // [[:], fai]

    // ASCAT reference files -- flat path collections (no meta wrapper), passed directly to ASCAT module
    allele_files = PREPARE_REFERENCE_FILES.out.allele_files  // [path, ...]  -- per-chromosome allele files
    loci_files   = PREPARE_REFERENCE_FILES.out.loci_files    // [path, ...]  -- per-chromosome loci files
    gc_file      = PREPARE_REFERENCE_FILES.out.gc_file       // [path, ...]  -- GC correction ([] if skipped)
    rt_file      = PREPARE_REFERENCE_FILES.out.rt_file       // [path, ...]  -- RT correction ([] if skipped)

    //
    // MODULE: FIBERTOOLSRS_PREDICTM6A
    //
    // predict m6a in unaligned bam

    if (!params.skip_fiber) {
        // Fiber-seq processing runs per-replicate on the unaligned BAMs.
        // Each replicate is processed independently; replicates are merged after alignment.
        if (!params.skip_normalfiber){
            // Process all samples (including normals) for fiber-seq
            ubams = ch_samplesheet
        }
        else {
            // Skip fiber-seq processing for normal samples; set aside normals to re-join later
            ch_samplesheet
            .branch { meta, _bams ->
                normal: meta.type == "normal"
                tumor: meta.type == "tumor"
                }
            .set { ch_ubams_normal_branching }
            // ch_ubams_normal_branching.normal: [meta, bam]  -- normal samples (held out)
            // ch_ubams_normal_branching.tumor:  [meta, bam]  -- tumor samples only

            normal_bams = ch_ubams_normal_branching.normal
            ubams = ch_ubams_normal_branching.tumor
        }
            // Branch by sequencing platform: PacBio needs m6A prediction, ONT does not
            ubams
            .branch{ meta, _bams ->
                pacBio: meta.platform == "pb"
                ont: meta.platform == "ont"
            }
            .set{ch_ubams_pacbio_ont_branching}
        // ch_ubams_pacbio_ont_branching.pacBio: [meta, bam]  -- PacBio samples
        // ch_ubams_pacbio_ont_branching.ont:    [meta, bam]  -- ONT samples (skip m6A)

        pacbio_bams = ch_ubams_pacbio_ont_branching.pacBio
        // Branch PacBio samples: only those with kinetics tags can have m6A predicted
        pacbio_bams
            .branch{meta, _bams ->
                kinetics: meta.kinetics == "true"
                noKinetics: meta.kinetics == "false"
            }
            .set{pacbio_bams}
        // pacbio_bams.kinetics:   [meta, bam]  -- PacBio with kinetics (mm/ml tags); m6A predictable
        // pacbio_bams.noKinetics: [meta, bam]  -- PacBio without kinetics; skip PREDICTM6A

        if (!params.skip_m6a) {
            //
            // MODULE: FIBERTOOLSRS_PREDICTM6A (label: process_high)
            // Input:  [meta, bam]  -- PacBio BAM with kinetics tags
            // Output: .bam -- [meta, bam]  -- BAM with m6A (MM/ML) tags added
            //
            FIBERTOOLSRS_PREDICTM6A (
                pacbio_bams.kinetics
            )
            // Merge PacBio with and without kinetics: both now have (or skip) m6A tags
            pacbio_bams.noKinetics
                .mix(FIBERTOOLSRS_PREDICTM6A.out.bam)
                .set{predicted_bams}
        }
        else {
            pacbio_bams.noKinetics
                .mix(pacbio_bams.kinetics)
                .set{predicted_bams}
        }
        // predicted_bams: [meta, bam]  -- all PacBio samples (m6A tags present where applicable)

        // Re-merge ONT and PacBio before fiber-seq branching
        ch_ubams_pacbio_ont_branching.ont
            .mix(predicted_bams)
            .set{fiber_branch}
        // fiber_branch (pre-split): [meta, bam]  -- all samples (ONT + PacBio, with m6A if applicable)

        // Branch on fiber-seq flag: only fiber-seq samples get nucleosome/FIRE calling
        fiber_branch
            .branch{ meta, _bams ->
                fiber: meta.fiber == "y"
                nonFiber: meta.fiber == "n"
            }
            .set{fiber_branch}
        // fiber_branch.fiber:    [meta, bam]  -- fiber-seq samples → nucleosome + FIRE calling
        // fiber_branch.nonFiber: [meta, bam]  -- non-fiber samples → passed through unchanged

        //
        // MODULE: FIBERTOOLSRS_NUCLEOSOMES (label: process_high)
        // Input:  [meta, bam]  -- fiber-seq BAM (with m6A tags for PacBio)
        // Output: .bam -- [meta, bam]  -- BAM with nucleosome footprint tags added
        //

        FIBERTOOLSRS_NUCLEOSOMES (
            fiber_branch.fiber
        )

        //
        // MODULE: FIBERTOOLSRS_FIRE (label: process_high)
        // Input:  FIBERTOOLSRS_NUCLEOSOMES.out.bam -- [meta, bam]  -- BAM with nucleosome tags
        // Output: .bam -- [meta, bam]  -- BAM with FIRE (Fiber-seq Inferred Regulatory Elements) tags
        //

        FIBERTOOLSRS_FIRE (
            FIBERTOOLSRS_NUCLEOSOMES.out.bam
        )

        if (!params.skip_normalfiber){
            // Re-merge fiber and non-fiber samples after FIRE annotation
            fiber_branch.nonFiber
            .mix(FIBERTOOLSRS_FIRE.out.bam)
            .set{ch_ubams}
        }
        else {
            // Re-merge fiber, non-fiber, and held-out normal samples
            fiber_branch.nonFiber
            .mix(normal_bams)
            .mix(FIBERTOOLSRS_FIRE.out.bam)
            .set{ch_ubams}
        }
        // ch_ubams (updated): [meta_with_replicate, bam]  -- all samples per-replicate;
        //   fiber-seq samples now carry nucleosome + FIRE tags; m6A tags for PacBio fiber-seq

        if(!params.skip_qc) {
            //
            // MODULE: FIBERTOOLSRS_QC (label: process_medium)
            // Input:  FIBERTOOLSRS_FIRE.out.bam -- [meta, bam]  -- annotated fiber-seq BAM
            // Output: QC reports for fiber-seq signal (written to outdir)
            //

            FIBERTOOLSRS_QC (
                FIBERTOOLSRS_FIRE.out.bam
            )
        }
    }
    //
    // MODULE: MINIMAP2_ALIGN (label: process_high)
    // Runs once per replicate.  The @RG header line encodes sample, type, and replicate so
    // that reads remain distinguishable after the per-replicate BAMs are merged below.
    // Input:  [meta_with_replicate, bam]  -- unaligned BAM per replicate
    //         ch_fasta     -- [[:], fasta]
    //         sort_bam=true, cigar_paf_format='bai', cigar_bam='', split_prefix=''
    // Output: .bam   -- [meta_with_replicate, bam]  -- coordinate-sorted aligned BAM
    //         .index -- [meta_with_replicate, bai]  -- BAM index
    //

    MINIMAP2_ALIGN (
        ch_ubams,
        ch_fasta,
        true,
        'bai',
        "",
        ""
    )

    // Join per-replicate BAM with its index, drop the replicate field, then group replicates
    // belonging to the same sample.  Single-replicate samples skip SAMTOOLS_MERGE.
    MINIMAP2_ALIGN.out.bam
        .join(MINIMAP2_ALIGN.out.index)
        .map { meta, bam, bai ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'type',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics',
                            'n_replicates')
            // groupKey tells groupTuple() how many items to expect for this group,
            // allowing it to release each sample as soon as all its replicates arrive
            // rather than waiting for all samples globally.
            return [groupKey(new_meta, new_meta.n_replicates), bam, bai]
        }
        .groupTuple()
        .map { meta, bams, bais ->
            [meta, bams.flatten(), bais.flatten()]
        }
        .branch { _meta, bams, _bais ->
            single:   bams.size() == 1
            multiple: bams.size() > 1
        }
        .set { ch_aligned_split }
    // ch_aligned_split.single:   [meta, [bam], [bai]]  -- one replicate; pass through
    // ch_aligned_split.multiple: [meta, [bam...], [bai...]]  -- merge needed

    // Single-replicate: unwrap lists to scalar paths
    ch_aligned_split.single
        .map { meta, bams, bais -> [meta, bams[0], bais[0]] }
        .set { ch_single_indexed }

    //
    // MODULE: SAMTOOLS_MERGE (label: process_low)
    // Merges per-replicate coordinate-sorted BAMs.  Because each was aligned with a unique
    // @RG line (ID = {sample}_{type}_rep{N}), the merged BAM retains full replicate identity.
    // Input:  [meta, [bam...], [bai...]]  -- grouped replicate BAMs + indices
    // Output: .bam -- [meta, bam]  -- merged BAM
    //
    SAMTOOLS_MERGE(
        ch_aligned_split.multiple,
        [[],[],[],[]]
    )

    // Index the merged BAM to produce a BAI (SAMTOOLS_MERGE does not create BAI inline)
    SAMTOOLS_INDEX_MERGE(SAMTOOLS_MERGE.out.bam)

    // Combine single-replicate and merged paths into a unified [meta, bam, bai] channel
    ch_single_indexed
        .mix(
            SAMTOOLS_MERGE.out.bam
                .join(SAMTOOLS_INDEX_MERGE.out.bai)
        )
        .set { ch_index_minimap }
    // ch_index_minimap: [meta, bam, bai]  -- one aligned BAM + index per sample (all replicates merged)

    // Convenience channel used by CRAMINO_POST and other modules that only need the BAM
    ch_index_minimap
        .map { meta, bam, _bai -> [meta, bam] }
        .set { ch_minimap_bam }
    // ch_minimap_bam: [meta, bam]  -- post-alignment BAM (replicates merged)

    //
    // MODULE: MODKIT_PILEUP
    //

    if (!params.skip_modkit) {
        MODKIT_PILEUP(ch_index_minimap, ch_fasta, ch_fai, [[:],[]])
    }

    ch_index_minimap
        .branch { meta, _bams, _bais ->
                paired: meta.paired_data       // meta.paired_data is the normal sample ID for tumors, or the tumor ID for normals
                tumor_only: !meta.paired_data  // meta.paired_data is null/false for tumor-only samples
        }
        .set { branched_minimap }

    // branched_minimap.paired:     [meta, bam, bai]  -- tumor AND normal samples flow together here;
    //                                                    each item is a single sample, joined downstream
    // branched_minimap.tumor_only: [meta, bam, bai]  -- tumor-only samples (no matched normal)

    // SUBWORKFLOW: TUMORONLY_SMALLVAR
    // Input:  branched_minimap.tumor_only -- [meta, bam, bai]
    // Output: .somatic_vcf  -- [meta, vcf, tbi]  -- somatic SNVs/indels
    //         .germline_vcf -- [meta, vcf, tbi]  -- germline SNVs/indels (ClairS-TO germline output)
    TUMORONLY_SMALLVAR(
        branched_minimap.tumor_only,
        ch_fasta,
        ch_fai,
        clairsto_pon_channel,
        ds_pon_channel
    )

    branched_minimap.paired
        .set{paired_ch}

    // Split paired samples into tumor and normal streams for joining
    paired_ch
        .branch { meta, _bams, _bais ->
                normal: meta.type == "normal"
                tumor:  meta.type == "tumor"
        }
        .set{branched_paired_ch}
    // branched_paired_ch.normal: [meta, bam, bai]  -- normal samples (meta.type == "normal")
    // branched_paired_ch.tumor:  [meta, bam, bai]  -- tumor samples  (meta.type == "tumor")

    // Strip 'type' field from normal meta before joining, so the key is just sample ID
     branched_paired_ch.normal
        .map{ meta, bam, bai ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, bam, bai]
        }
        .set{paired_normal_bams}
    // paired_normal_bams: [meta (no type), normal_bam, normal_bai]

    // Join tumor and normal BAMs into a single channel for somatic variant calling
    // Join key is meta (with 'type' stripped), so tumor meta.id must equal normal meta.id
    branched_paired_ch.tumor
        .map{ meta, bam, bai ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, bam, bai]
        }
        .join(paired_normal_bams)
        .set { somatic_smallvar_input }
    // somatic_smallvar_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

    // SUBWORKFLOW: PAIRED_SMALLVAR_SOMATIC
    // Input:  somatic_smallvar_input -- [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
    // Output: .somatic_vcf -- [meta, vcf, tbi]  -- somatic SNVs/indels (ClairS and/or DeepSomatic consensus)
    PAIRED_SMALLVAR_SOMATIC (
        somatic_smallvar_input,
        ch_fasta,
        ch_fai,
        ds_pon_channel
    )

    // SUBWORKFLOW: PAIRED_SMALLVAR_GERMLINE
    // Input:  branched_paired_ch.normal -- [meta, bam, bai]  -- normal sample BAMs only
    //         downloaded_clair3_models  -- [meta(id=model_name), model_dir]
    // Output: .germline_vcf -- [meta, vcf, tbi]  -- germline SNVs/indels (Clair3 and/or DeepVariant consensus)
    PAIRED_SMALLVAR_GERMLINE (
        branched_paired_ch.normal,
        ch_fasta,
        ch_fai,
        downloaded_clair3_models
    )

    // Merge germline VCFs from paired and tumor-only paths into a single channel
    PAIRED_SMALLVAR_GERMLINE.out.germline_vcf
        .mix(TUMORONLY_SMALLVAR.out.germline_vcf)
        .set{ch_germline_vcf}
    // ch_germline_vcf: [meta, vcf, tbi]  -- germline variants for all samples (paired + tumor-only)

    // Merge somatic VCFs from tumor-only and paired T/N paths into a single channel
    TUMORONLY_SMALLVAR.out.somatic_vcf
        .mix(PAIRED_SMALLVAR_SOMATIC.out.somatic_vcf)
        .set{ch_somatic_vcf}
    // ch_somatic_vcf: [meta, vcf, tbi]  -- somatic variants for all samples

    // SUBWORKFLOW: PHASING_HAPLOTYPING
    // Input:  ch_index_minimap -- [meta, bam, bai]  -- all aligned BAMs (tumor + normal + tumor-only)
    //         ch_germline_vcf  -- [meta, vcf, tbi]  -- germline variants (used to phase reads)
    //         ch_somatic_vcf   -- [meta, vcf, tbi]  -- somatic variants (get phasing transferred)
    //         ch_fasta / ch_fai
    // Output: .phased_germline_vcf    -- [meta, vcf, tbi]  -- phased germline VCF
    //         .phased_somatic_vcf     -- [meta, vcf, tbi]  -- phased somatic VCF
    //         .tumor_normal_hapbams_ch -- [meta, bam, bai] -- haplotagged BAMs (all samples)
    PHASING_HAPLOTYPING (
        ch_index_minimap,
        ch_germline_vcf,
        ch_somatic_vcf,
        ch_fasta,
        ch_fai
    )

    // Prepare phased VCFs for VEP: add empty 'extra' list required by ENSEMBLVEP_VEP
    PHASING_HAPLOTYPING.out.phased_somatic_vcf
        .map { meta, vcf, _tbi ->
            def extra = []
            return [meta, vcf, extra]
        }
        .set { somatic_vep }
    // somatic_vep: [meta, vcf, []]  -- phased somatic VCF ready for VEP annotation

    PHASING_HAPLOTYPING.out.phased_germline_vcf
        .map { meta, vcf, _tbi ->
            def extra = []
            return [meta, vcf, extra]
        }
        .set { germline_vep }
    // germline_vep: [meta, vcf, []]  -- phased germline VCF ready for VEP annotation


    whatshap_stats_txt = channel.empty()

    if (!params.skip_qc && !params.skip_whatshapstats) {

        // Drop the empty 'extra' element added for VEP input
        germline_vep
            .map { meta, vcf, _extra ->
                return [meta, vcf] }
            .set { ch_whatshap_stats }
        // ch_whatshap_stats: [meta, vcf]  -- phased germline VCF for phasing QC

        //
        // MODULE: WHATSHAP_STATS (label: process_single)
        // Input:  [meta, vcf]   -- phased VCF (germline)
        //         gtf=true, sample=true, chr_lengths=false
        // Output: .tsv -- [meta, tsv]  -- per-chromosome phasing statistics
        //

        WHATSHAP_STATS (
            ch_whatshap_stats,
            true,
            true,
            false
        )

        whatshap_stats_txt = WHATSHAP_STATS.out.tsv

    }

    if (!params.skip_vep) {

        //
        // MODULE: GERMLINE_VEP (ENSEMBLVEP_VEP alias; label: process_medium)
        // Input:  germline_vep -- [meta, vcf, []]  -- phased germline VCF
        //         vep_cache    -- [[:], cache_dir]
        //         ch_fasta     -- [[:], fasta]
        // Output: annotated germline VCF with consequence predictions
        //
        if (params.vep_custom != null) {
            vep_custom = file(params.vep_custom)
        } else {
            vep_custom = []
        }
        if (params.vep_custom_tbi != null) {
            vep_custom_tbi = file(params.vep_custom_tbi)
        } else {
            vep_custom_tbi = []
        }
        GERMLINE_VEP (
            germline_vep,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache,
            ch_fasta,
            [],
            vep_custom,
            vep_custom_tbi
        )

        //
        // MODULE: SOMATIC_VEP (ENSEMBLVEP_VEP alias; label: process_medium)
        // Input:  somatic_vep -- [meta, vcf, []]  -- phased somatic VCF
        //         vep_cache   -- [[:], cache_dir]
        //         ch_fasta    -- [[:], fasta]
        // Output: annotated somatic VCF with consequence predictions
        //

        SOMATIC_VEP (
            somatic_vep,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache,
            ch_fasta,
            [],
            vep_custom,
            vep_custom_tbi
        )
    }

    // Build SEVERUS input by combining tumor-only and T/N paired samples with phased germline VCFs
    // Tumor-only samples get empty lists for normal BAM/BAI (SEVERUS runs in tumor-only mode)
    branched_minimap.tumor_only
        .map{ meta, bam, bai ->
            def new_meta = meta.subMap('id',
                            'paired_data',
                            'platform',
                            'sex',
                            'fiber',
                            'clair3_model',
                            'clairS_model',
                            'clairSTO_model',
                            'kinetics')
            return[new_meta, bam, bai]
        }
        .map{meta, tumor_bam, tumor_bai->
            def normal_bam = []
            def normal_bai = []
            return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
        }
        // Mix with paired T/N input (which already has normal BAM/BAI from somatic_smallvar_input)
        .mix(somatic_smallvar_input)
        // Attach phased germline VCF (used by SEVERUS for phased SV calling)
        .join(PHASING_HAPLOTYPING.out.phased_germline_vcf)
        .set{severus_input}
    // severus_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_germline_vcf, phased_germline_tbi]
    //   normal_bam/bai are empty lists [] for tumor-only samples

    //
    // MODULE: SEVERUS (label: process_high)
    // Input:  severus_input -- [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, vcf, tbi]
    //         [[:], bed_file, pon_file]  -- optional target BED and panel-of-normals for SV filtering
    // Output: .all_vcf -- [meta, vcf]  -- all somatic SVs (sniffles2 format)
    //

    SEVERUS (
        severus_input,
        [[:], params.bed_file, params.pon_file]
    )

    ch_versions = ch_versions.mix(SEVERUS.out.versions)

    SEVERUS.out.all_vcf
        .map { meta, vcf ->
            def extra = []
            return [meta, vcf, extra]
        }
        .set { sv_vep }
    // sv_vep: [meta, severus_all_vcf, []]  -- all SVs ready for VEP annotation

    if(!params.skip_vep) {
        //
        // MODULE: SV_VEP (ENSEMBLVEP_VEP alias; label: process_medium)
        // Input:  sv_vep -- [meta, vcf, []]  -- SEVERUS SV VCF
        // Output: annotated SV VCF with consequence predictions
        //
        SV_VEP (
            sv_vep,
            params.vep_genome,
            params.vep_species,
            params.vep_cache_version,
            vep_cache,
            ch_fasta,
            [],
            vep_custom,
            vep_custom_tbi
        )
    }


    ch_nanoplot_post_txt = channel.empty()


    if (!params.skip_qc && !params.skip_cramino) {

        //
        // MODULE: CRAMINO_POST (label: process_medium)
        // Input:  ch_minimap_bam -- [meta, bam]  -- post-alignment coordinate-sorted BAM
        // Output: .arrow -- [meta, arrow_file]  -- alignment statistics in feather format
        //

        CRAMINO_POST ( ch_minimap_bam )

        if (!params.skip_nanoplot) {

            //
            // MODULE: NANOPLOT_POST (label: process_medium)
            // Input:  CRAMINO_POST.out.arrow -- [meta, arrow_file]
            // Output: HTML/txt QC reports (post-alignment)
            //

            NANOPLOT_POST(CRAMINO_POST.out.arrow)

        }


    }

    //
    // Module: MOSDEPTH
    //

    ch_mosdepth_global = channel.empty()
    ch_mosdepth_summary = channel.empty()

    if (!params.skip_qc && !params.skip_mosdepth) {

        // MOSDEPTH requires a BED file argument; pass [] to compute genome-wide depth
        ch_index_minimap
            .map { meta, bam, bai -> [meta, bam, bai, []] }
            .set { ch_mosdepth_in }
        // ch_mosdepth_in: [meta, bam, bai, []]  -- [] is the optional BED (empty = genome-wide)

        //
        // MODULE: MOSDEPTH (label: process_medium)
        // Input:  [meta, bam, bai, bed]  -- bed is [] for genome-wide coverage
        //         ch_fasta -- [[:], fasta]  -- used for CRAM decoding (if applicable)
        // Output: .global_txt  -- [meta, txt]  -- global depth summary
        //         .summary_txt -- [meta, txt]  -- per-contig depth summary
        //
        MOSDEPTH (
            ch_mosdepth_in,
            ch_fasta
        )

        ch_mosdepth_global = MOSDEPTH.out.global_txt
        ch_mosdepth_summary = MOSDEPTH.out.summary_txt
    }

    //
    // SUBWORKFLOW: BAM_STATS_SAMTOOLS (nf-core subworkflow)
    // Input:  [meta, bam, bai]  -- aligned BAM with index
    //         ch_fasta          -- [[:], fasta]
    // Output: .stats    -- [meta, txt]  -- samtools stats output
    //         .flagstat -- [meta, txt]  -- samtools flagstat output
    //         .idxstats -- [meta, txt]  -- samtools idxstats output
    //
    ch_bam_stats = channel.empty()
    ch_bam_flagstat = channel.empty()
    ch_bam_idxstats = channel.empty()

    if (!params.skip_qc && !params.skip_bamstats ) {

        BAM_STATS_SAMTOOLS (
            ch_index_minimap, // [meta, bam, bai]
            ch_fasta
        )

        ch_bam_stats = BAM_STATS_SAMTOOLS.out.stats
        ch_bam_flagstat = BAM_STATS_SAMTOOLS.out.flagstat
        ch_bam_idxstats = BAM_STATS_SAMTOOLS.out.idxstats
    }

    //
    // MODULE: ASCAT (label: process_high)
    // Input:  [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]  -- NOTE: normal before tumor (ASCAT convention)
    //         allele_files, loci_files, gc_file, rt_file  -- ASCAT reference files
    // Output: .png plots, .segments, .purity_ploidy  -- copy number results
    //

    if (!params.skip_ascat) {
        // ASCAT expects [normal, tumor] order; rearrange from severus_input [tumor, normal] order
        severus_input
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, _vcf, _tbi ->
                return [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]
            }
            .set { ascat_ch }
        // ascat_ch: [meta, normal_bam, normal_bai, tumor_bam, tumor_bai]

        ASCAT (
            ascat_ch,
            params.genome_name,
            allele_files,
            loci_files,
            [],
            [],
            gc_file,
            rt_file
        )

        ch_versions = ch_versions.mix(ASCAT.out.versions)
    }

    //
    // SUBWORKFLOWS: TUMORONLY_SAVANA / PAIRED_SAVANA (SAVANA SV + copy-number calling)
    // Tumor-only runs the combined `savana to`; matched tumor/normal runs `savana run` +
    // `savana classify` + `savana cna` as separate steps (mirrors Severus/small-variant split).
    // CN calls are auto-published per-process (conf/modules.config); somatic_vcf is fed into
    // SV_VEP below, alongside Severus's SVs.
    //

    savana_somatic_vcf = channel.empty()

    if (!params.skip_savana) {
        // SAVANA reads the HP (haplotype) tag per read and its README recommends phased BAMs,
        // so build its input from PHASING_HAPLOTYPING's haplotagged BAMs rather than the
        // unphased ones severus_input carries.
        def savana_meta_keys = ['id', 'paired_data', 'platform', 'sex', 'fiber',
            'clair3_model', 'clairS_model', 'clairSTO_model', 'kinetics']

        PHASING_HAPLOTYPING.out.tumor_normal_hapbams_ch
            .branch { meta, _bam, _bai ->
                tumor_only:    !meta.paired_data
                paired_tumor:  meta.paired_data && meta.type == 'tumor'
                paired_normal: meta.paired_data && meta.type == 'normal'
            }
            .set { branched_savana_hapbams }

        branched_savana_hapbams.tumor_only
            .map { meta, bam, bai ->
                def normal_bam = []
                def normal_bai = []
                return [meta.subMap(savana_meta_keys), bam, bai, normal_bam, normal_bai]
            }
            .set { savana_tumoronly_hapbams }

        branched_savana_hapbams.paired_tumor
            .map { meta, bam, bai -> return [meta.subMap(savana_meta_keys), bam, bai] }
            .set { savana_paired_tumor_hapbams }

        branched_savana_hapbams.paired_normal
            .map { meta, bam, bai -> return [meta.subMap(savana_meta_keys), bam, bai] }
            .set { savana_paired_normal_hapbams }

        savana_paired_tumor_hapbams
            .join(savana_paired_normal_hapbams)
            .set { savana_paired_hapbams }
        // savana_paired_hapbams: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]

        savana_tumoronly_hapbams
            .mix(savana_paired_hapbams)
            .join(PHASING_HAPLOTYPING.out.phased_germline_vcf)
            .set { savana_input }
        // savana_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_germline_vcf, phased_germline_tbi]
        //   normal_bam/bai are [] for tumor-only samples

        savana_input
            .branch { meta, _tumor_bam, _tumor_bai, normal_bam, _normal_bai, _phased_vcf, _phased_tbi ->
                tumor_only: !normal_bam
                paired: normal_bam
            }
            .set { branched_savana_input }
        // branched_savana_input.tumor_only / .paired: same 7-tuple shape as savana_input

        branched_savana_input.tumor_only
            .map { meta, tumor_bam, tumor_bai, _normal_bam, _normal_bai, phased_vcf, phased_tbi ->
                return [meta, tumor_bam, tumor_bai, phased_vcf, phased_tbi]
            }
            .set { tumoronly_savana_input }
        // tumoronly_savana_input: [meta, tumor_bam, tumor_bai, phased_vcf, phased_tbi]

        ch_savana_contigs = channel.value([[:], params.savana_contigs])
        // Tumor-only has no matched germline control, so allele counting uses the bundled 1000g
        // population SNP set instead of a (nonexistent) germline VCF -- see TUMORONLY_SAVANA.
        ch_savana_g1000_vcf = channel.value([[:], params.savana_g1000_vcf])

        TUMORONLY_SAVANA (
            tumoronly_savana_input,
            ch_fasta,
            ch_fai,
            ch_savana_contigs,
            ch_savana_g1000_vcf
        )

        PAIRED_SAVANA (
            branched_savana_input.paired,
            ch_fasta,
            ch_fai,
            ch_savana_contigs
        )

        TUMORONLY_SAVANA.out.somatic_vcf
            .mix(PAIRED_SAVANA.out.somatic_vcf)
            .set { savana_somatic_vcf }
        // savana_somatic_vcf: [meta, vcf]

        if (!params.skip_vep) {
            //
            // MODULE: VEP_SAVANA (ENSEMBLVEP_VEP alias; label: process_medium)
            // Input:  savana_somatic_vcf -- [meta, vcf, []]  -- SAVANA classified somatic SV VCF
            // Output: annotated SV VCF with consequence predictions
            //
            savana_somatic_vcf
                .map { meta, vcf ->
                    def extra = []
                    return [meta, vcf, extra]
                }
                .set { savana_vep }
            // savana_vep: [meta, savana_somatic_vcf, []]  -- SAVANA SVs ready for VEP annotation

            VEP_SAVANA (
                savana_vep,
                params.vep_genome,
                params.vep_species,
                params.vep_cache_version,
                vep_cache,
                ch_fasta,
                [],
                vep_custom,
                vep_custom_tbi
            )
        }
    }

    //
    // MODULE: WAKHAN (label: process_medium)
    // Haplotype-aware genome assembly and variant phasing visualisation
    // Input:  [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_germline_vcf, severus_all_vcf]
    //         ch_fasta          -- [[:], fasta]
    //         centromere_bed    -- BED file of centromere coordinates (for assembly anchoring)
    // Output: WAKHAN assembly reports (written to outdir)
    //

    if (!params.skip_wakhan) {

        // Attach SEVERUS SV VCF to the severus_input channel (dropping the phased TBI)
        severus_input
            .join(SEVERUS.out.all_vcf)
            .map { meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, _phased_tbi, all_vcf ->
                return [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_vcf, all_vcf]
            }
            .set { wakhan_input }
        // wakhan_input: [meta, tumor_bam, tumor_bai, normal_bam, normal_bai, phased_germline_vcf, severus_all_vcf]
        //   normal_bam/bai are [] for tumor-only samples

        WAKHAN (
            wakhan_input,
            ch_fasta,
            file(params.centromere_bed)
        )
    }

    //
    // Collate software versions from two sources:
    //   1. ch_versions (classic path): version YAML files emitted by modules
    //   2. channel.topic("versions") (topic channel path): version tuples [process, tool, version]
    //      emitted directly by modules that use the topic-channel pattern
    //
    def topic_versions = channel.topic("versions")
        .distinct()  // deduplicate identical version entries across samples
        .branch { entry ->
            versions_file:  entry instanceof Path   // classic YAML file path
            versions_tuple: true                    // [process, tool, version] tuple
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            // Strip workflow prefix (everything before the last ':') from process name
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)  // group tool versions by process name
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
    // topic_versions_string: formatted YAML-like string per process, ready to write

    // Merge both version sources and write to versions YAML (consumed by MultiQC)
    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'lrsomatic_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
    // ch_collated_versions: path  -- merged software versions YAML for MultiQC


    //
    // MODULE: MULTIQC (label: process_single)
    // Aggregates QC reports from all modules into a single HTML report
    // Input:  [[id:'multiqc'], [qc_files...], [config_files...], [logo], [], []]
    // Output: .report -- [meta, html]  -- MultiQC HTML report
    //
    summary_params = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    // Collect QC outputs from all optional modules
    // .collect{it -> it[1]} extracts the file from [meta, file] tuples; ifEmpty([]) handles skipped modules
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_stats.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_flagstat.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_bam_idxstats.collect{it -> it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(ch_mosdepth_global.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_mosdepth_summary.collect{it -> it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(ch_nanoplot_pre_txt.collect{it -> it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_nanoplot_post_txt.collect{it -> it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(whatshap_stats_txt.collect{it -> it[1]}.ifEmpty([]))

    // Build the final MULTIQC input tuple: all QC files + config files + logo
    MULTIQC (
        ch_multiqc_files
            .collect()
            .map { files ->
                def multiqc_config_files = [file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)]
                if (params.multiqc_config) {
                    multiqc_config_files += [file(params.multiqc_config, checkIfExists: true)]
                }
                def multiqc_logo_file = params.multiqc_logo ? [file(params.multiqc_logo, checkIfExists: true)] : []
                // MULTIQC input: [meta, [qc_files], [config_files], [logo], [], []]
                [[id: 'multiqc'], files, multiqc_config_files, multiqc_logo_file, [], []]
            }
    )

    emit:
        multiqc_report = MULTIQC.out.report.map { _meta, report -> report } // channel: /path/to/multiqc_report.html
        versions       = ch_versions                 // channel: [ path(versions.yml) ]



}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
