include { SIGPROFILER_INSTALL } from '../../modules/local/sigprofiler/install/main'

workflow PREPARE_SIGNATURES {

    take:
        genome              // str:  SigProfilerMatrixGenerator genome name (e.g. "GRCh38", "CHM13-T2T")
        genome_url          // str:  URL of <genome>.tar.gz, or null to let SigProfilerMatrixGenerator download from the AlexandrovLab FTP
        genome_dir          // path: existing SigProfilerMatrixGenerator volume (contains tsb/<genome>/), or null
        download_genome     // bool: if true, install the genome payload with SIGPROFILER_INSTALL instead of using genome_dir

    main:

        ch_versions = channel.empty()
        sigprofiler_volume = channel.empty()

        if (!genome) {
            error("No SigProfilerMatrixGenerator genome is defined for --genome ${params.genome}. Set --sigprofiler_genome (e.g. GRCh38 or CHM13-T2T) or use --skip_signatures.")
        }

        if (download_genome) {
            //
            // MODULE: SIGPROFILER_INSTALL (label: process_single, process_long)
            // Downloads and verifies the per-genome TSB payload (~3 GB) into a volume directory.
            // Published to ${params.outdir}/cache/ so later runs can pass it as --sigprofiler_genome_dir.
            //
            SIGPROFILER_INSTALL (
                genome,
                genome_url ?: ''
            )
            sigprofiler_volume = SIGPROFILER_INSTALL.out.volume
        }
        else {
            if (!genome_dir) {
                error("Mutational signature analysis needs the SigProfilerMatrixGenerator payload for ${genome}.\n" +
                      "Either pass --sigprofiler_genome_dir <volume> (a directory containing tsb/${genome}/) " +
                      "or add --download_sigprofiler_genome to install it into ${params.outdir}/cache/ on this run, " +
                      "or disable the step with --skip_signatures.")
            }
            def tsb_dir = file("${genome_dir}/tsb/${genome}", type: 'dir')
            if (!tsb_dir.exists() || !tsb_dir.isDirectory()) {
                error("Path provided with --sigprofiler_genome_dir is invalid.\nMake sure there is a directory named tsb/${genome} in ${genome_dir}.")
            }
            def n_chrom = tsb_dir.listFiles().count { f -> f.name.endsWith('.txt') }
            if (n_chrom < 24) {
                error("${tsb_dir} holds ${n_chrom} chromosome files; a complete SigProfilerMatrixGenerator install of ${genome} has 24 (1-22, X, Y).")
            }
            sigprofiler_volume = channel.value(file(genome_dir, type: 'dir', checkIfExists: true))
        }
        // sigprofiler_volume: path -- SigProfilerMatrixGenerator volume root (downloaded or validated local)

    emit:
        volume   = sigprofiler_volume  // path -- volume directory containing tsb/<genome>/
        versions = ch_versions
}
