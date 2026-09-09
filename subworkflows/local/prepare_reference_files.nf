//
// Prepare reference files (unzipping and adding index)
//

include { PIGZ_UNCOMPRESS as UNZIP_FASTA } from '../../modules/nf-core/pigz/uncompress/main'
include { SAMTOOLS_FAIDX                 } from '../../modules/nf-core/samtools/faidx/main'
include { UNZIP as UNZIP_ALLELES         } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_GC              } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_LOCI            } from '../../modules/nf-core/unzip/main'
include { UNZIP as UNZIP_RT              } from '../../modules/nf-core/unzip/main'
include { UNTAR                          } from '../../modules/nf-core/untar/main'
include { WGET                           } from '../../modules/nf-core/wget/main'

workflow PREPARE_REFERENCE_FILES {
    take:
        fasta           // str: path to reference FASTA (may be .gz)
        ascat_alleles   // str: path to ASCAT allele files (directory or .zip), or null
        ascat_loci      // str: path to ASCAT loci files (directory or .zip), or null
        ascat_loci_gc   // str: path to ASCAT GC correction file (.zip or direct), or null
        ascat_loci_rt   // str: path to ASCAT RT correction file (.zip or direct), or null
        basecall_meta   // [meta, basecall_model_str, kinetics_str]  -- from METAEXTRACT per sample
        clair3_modelMap // Map<basecall_model_str, clair3_model_name>  -- used to resolve download URLs

    main:
        ch_versions = channel.empty()
        ch_prepared_fasta = channel.empty()
        allele_files = channel.empty()
        loci_files = channel.empty()
        gc_file = channel.empty()
        rt_file = channel.empty()

        // Decompress FASTA if gzipped; pass through as-is if already uncompressed
        if (fasta.endsWith('.gz')){
            //
            // MODULE: UNZIP_FASTA (PIGZ_UNCOMPRESS alias; label: process_medium)
            // Input:  [[:], fasta.gz]
            // Output: .file -- [[:], fasta]  -- decompressed FASTA
            //
            UNZIP_FASTA( [ [:], fasta ])

            ch_prepared_fasta = UNZIP_FASTA.out.file
            ch_versions = ch_versions.mix(UNZIP_FASTA.out.versions)
        } else {
            ch_prepared_fasta = channel.value([ [:], fasta ])
        }
        // ch_prepared_fasta: [[:], fasta_path]  -- empty meta; uncompressed FASTA

        // Build Clair3 model download URLs from basecall metadata
        // Priority: explicit meta.clair3_model param > auto-detected from BAM header via modelMap
        // PacBio models from HKU mirror; ONT models from Oxford Nanopore CDN
        basecall_meta.map { meta, basecall_model_meta, _kinetics_meta ->
            def model = (!meta.clair3_model || meta.clair3_model.toString().trim() in ['', '[]']) ? clair3_modelMap.get(basecall_model_meta) : meta.clair3_model
            // Key the entry on the model that is actually downloaded. Keying on the header-derived
            // name instead made a sample with an explicit clair3_model produce a second entry under
            // the wrong name: .unique() kept both, UNTAR extracted two different models into
            // directories with the same name, and the by-name combine in PAIRED_SMALLVAR_GERMLINE
            // ran Clair3 twice per normal BAM (once with the wrong model), with the downstream join
            // taking whichever finished first.
            def meta_new = [id: model]
            def download_prefix = ( basecall_model_meta == 'hifi_revio' ? "https://www.bio8.cs.hku.hk/clair3/clair3_models/" : "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3" )
            def url = "${download_prefix}/${model}.tar.gz"
            return [ meta_new, url ]
        }
        .unique()  // deduplicate: multiple samples with the same Clair3 model share one download
        .set{ clair3_model_urls }
        // clair3_model_urls: [meta(id=clair3_model_name), download_url_str]
        //   one item per unique Clair3 model needed across all samples

        //
        // MODULE: WGET (label: process_single)
        // Input:  [meta, url_str]  -- model name (id) + download URL
        // Output: .outfile -- [meta, tarball]  -- downloaded .tar.gz model archive
        //
        WGET ( clair3_model_urls )

        ch_versions = ch_versions.mix(WGET.out.versions)

        //
        // MODULE: UNTAR (label: process_single)
        // Input:  WGET.out.outfile -- [meta, tarball]
        // Output: .untar -- [meta, model_dir]  -- extracted Clair3 model directory
        //
        UNTAR (
            WGET.out.outfile
        )

        UNTAR.out.untar.set { downloaded_clair3_models }
        // downloaded_clair3_models: [meta(id=clair3_model_name), model_dir]

        //
        // MODULE: SAMTOOLS_FAIDX (label: process_single)
        // Input:  [[:], fasta, []]  -- empty meta + empty regions file (index full FASTA)
        //         false             -- do not write fai to stdout
        // Output: .fai -- [[:], fai_path]
        //
        SAMTOOLS_FAIDX (
            ch_prepared_fasta.map { meta, fa -> [meta, fa, []] },
            false
        )

        ch_prepared_fai = SAMTOOLS_FAIDX.out.fai
        // ch_prepared_fai: [[:], fai_path]  -- empty meta

        //
        // Prepare ASCAT reference files
        // Each file set can be provided as a .zip archive or a plain directory/file path
        // All ASCAT outputs are flat file collections (no meta tuple) for use with ASCAT module
        //
        if ( !params.skip_ascat ) {
            // Allele files: per-chromosome SNP allele frequency files (used for LogR/BAF calculation)
            if (!ascat_alleles) allele_files = channel.empty()
            else if (ascat_alleles.endsWith(".zip")) {
                // MODULE: UNZIP_ALLELES (UNZIP alias; label: process_single)
                // Input:  [meta(id=basename), [zip_file]]  -- collected zip
                // Output: .unzipped_archive -- [meta, dir]  -- extracted directory; flatMap lists individual files
                UNZIP_ALLELES(channel.fromPath(file(ascat_alleles)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                allele_files = UNZIP_ALLELES.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // allele_files: [path, path, ...]  -- all per-chromosome allele files collected
                ch_versions = ch_versions.mix(UNZIP_ALLELES.out.versions)
            } else allele_files = channel.fromPath(ascat_alleles).collect()

            // Loci files: per-chromosome SNP loci positions
            if (!ascat_loci) loci_files = channel.empty()
            else if (ascat_loci.endsWith(".zip")) {
                // MODULE: UNZIP_LOCI (UNZIP alias; label: process_single)
                UNZIP_LOCI(channel.fromPath(file(ascat_loci)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                loci_files = UNZIP_LOCI.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // loci_files: [path, path, ...]  -- all per-chromosome loci files collected
                ch_versions = ch_versions.mix(UNZIP_LOCI.out.versions)
            } else loci_files = channel.fromPath(ascat_loci).collect()

            // GC correction file: genome-wide GC content per locus (optional)
            if (!ascat_loci_gc) gc_file = channel.value([])
            else if ( ascat_loci_gc.endsWith(".zip") ) {
                // MODULE: UNZIP_GC (UNZIP alias; label: process_single)
                UNZIP_GC(channel.fromPath(file(ascat_loci_gc)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                gc_file = UNZIP_GC.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // gc_file: [path, ...]  -- GC correction file(s) collected
                ch_versions = ch_versions.mix(UNZIP_GC.out.versions)
            } else gc_file = channel.fromPath(ascat_loci_gc).collect()

            // Replication timing correction file: RT correction per locus (optional)
            if (!ascat_loci_rt) rt_file = channel.value([])
            else if (ascat_loci_rt.endsWith(".zip")) {
                // MODULE: UNZIP_RT (UNZIP alias; label: process_single)
                UNZIP_RT(channel.fromPath(file(ascat_loci_rt)).collect().map{ it -> [ [ id:it[0].baseName ], it ] })
                rt_file = UNZIP_RT.out.unzipped_archive.flatMap { it -> it[1].listFiles() }.collect()
                // rt_file: [path, ...]  -- RT correction file(s) collected
                ch_versions = ch_versions.mix(UNZIP_RT.out.versions)
            } else rt_file = channel.fromPath(ascat_loci_rt).collect()
        }

    emit:
        prepped_fasta = ch_prepared_fasta  // [[:], fasta_path]  -- uncompressed reference FASTA
        prepped_fai   = ch_prepared_fai    // [[:], fai_path]    -- samtools FAI index

        // ASCAT reference files -- flat file collections (no meta tuple wrapper)
        // Each is a list of paths collected into a single channel value
        allele_files  // [path, ...]  -- per-chromosome allele frequency files
        loci_files    // [path, ...]  -- per-chromosome loci position files
        gc_file       // [path, ...]  -- GC correction file ([] if not provided)
        rt_file       // [path, ...]  -- replication timing correction file ([] if not provided)

        downloaded_clair3_models  // [meta(id=clair3_model_name), model_dir]

        versions = ch_versions
}
