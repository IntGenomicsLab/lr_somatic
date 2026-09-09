process SIGPROFILER_MATRIXGENERATOR {
    tag "$meta.id"
    label 'process_medium'

    // Conda is not supported: the image installs SigProfilerMatrixGenerator from the fork
    // that adds the CHM13-T2T genome (SigProfilerSuite/SigProfilerMatrixGenerator#250) and
    // SigProfilerAssignment from the fork that adds CHM13-T2T COSMIC signatures. Return to
    // the bioconda/biocontainers releases once both are merged and released upstream.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/sigprofiler-sif:1.3.6-chm13-28a9ce8'
        : 'ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-28a9ce8'}"

    input:
    tuple val(meta), path(vcf)          // somatic small-variant VCF (plain or bgzipped)
    // Staged under a fixed name: the script builds its own local `volume/`, and the directory
    // produced by SIGPROFILER_INSTALL (and the published cache) is itself called `volume`.
    tuple val(meta2), path(volume, stageAs: 'genome_volume')      // SigProfilerMatrixGenerator volume containing tsb/<genome>/
    val(genome)                         // SigProfilerMatrixGenerator genome name, e.g. GRCh38 or CHM13-T2T

    output:
    tuple val(meta), path("output/SBS/*.SBS96.all")  , emit: sbs96
    tuple val(meta), path("output/DBS/*.DBS78.all")  , emit: dbs78    , optional: true
    tuple val(meta), path("output/ID/*.ID83.all")    , emit: id83     , optional: true
    tuple val(meta), path("output")                  , emit: output_dir
    tuple val(meta), path("output/logs/*")           , emit: logs     , optional: true
    tuple val("${task.process}"), val('sigprofilermatrixgenerator'), eval("python -c 'import importlib.metadata as m; print(m.version(\"SigProfilerMatrixGenerator\"))'"), topic: versions, emit: versions_sigprofilermatrixgenerator

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SIGPROFILER_MATRIXGENERATOR does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // SigProfilerMatrixGenerator names the sample after the input file name up to the first '.',
    // and dispatches on the extension of the (single) file in the input directory.
    def sample = prefix.replaceAll(/[^A-Za-z0-9_-]/, '_')
    def decompress = vcf.name.endsWith('.gz') ? "gzip -cd ${vcf} > input/${sample}.vcf" : "cp ${vcf} input/${sample}.vcf"
    """
    # Stage only a symlink to the genome payload in a local volume so the tool never writes
    # into the shared reference directory.
    mkdir -p input volume/tsb plot_templates
    ln -s "\$(readlink -f ${volume})/tsb/${genome}" volume/tsb/${genome}
    # sigProfilerPlotting caches plot templates in its (read-only) package directory unless told otherwise
    export SIGPROFILERPLOTTING_VOLUME="\$PWD/plot_templates"

    ${decompress}

    SigProfilerMatrixGenerator matrix_generator \\
        ${prefix} \\
        ${genome} \\
        input/ \\
        --volume volume \\
        ${args}

    # SigProfilerMatrixGenerator writes everything (matrices, plots, logs, sorted VCF) under input/output
    mv input/output output
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p output/SBS output/DBS output/ID output/logs
    touch output/SBS/${prefix}.SBS96.all output/DBS/${prefix}.DBS78.all output/ID/${prefix}.ID83.all output/logs/${prefix}.out
    """
}
