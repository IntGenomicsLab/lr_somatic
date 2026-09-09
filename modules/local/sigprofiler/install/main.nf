process SIGPROFILER_INSTALL {
    tag "$genome"
    label 'process_single'
    label 'process_long'

    // Conda is not supported: the image installs SigProfilerMatrixGenerator from the fork
    // that adds the CHM13-T2T genome (SigProfilerSuite/SigProfilerMatrixGenerator#250) and
    // SigProfilerAssignment from the fork that adds CHM13-T2T COSMIC signatures. Return to
    // the bioconda/biocontainers releases once both are merged and released upstream.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/sigprofiler-sif:1.3.6-chm13-28a9ce8'
        : 'ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-28a9ce8'}"

    input:
    val(genome)      // SigProfilerMatrixGenerator genome name, e.g. GRCh38 or CHM13-T2T
    val(url)         // URL of <genome>.tar.gz to install from; empty/null = AlexandrovLab FTP

    output:
    path("volume")                                                                                                           , emit: volume
    tuple val("${task.process}"), val('sigprofilermatrixgenerator'), eval("python -c 'import importlib.metadata as m; print(m.version(\"SigProfilerMatrixGenerator\"))'"), topic: versions, emit: versions_sigprofilermatrixgenerator

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SIGPROFILER_INSTALL does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args    = task.ext.args ?: ''
    def install = url
        ? """
    python -c 'import urllib.request; urllib.request.urlretrieve("${url}", "${genome}.tar.gz")'
    SigProfilerMatrixGenerator install ${genome} --local_genome . --volume volume ${args}
    rm -f ${genome}.tar.gz
    """
        : """
    SigProfilerMatrixGenerator install ${genome} --volume volume ${args}
    """
    """
    mkdir -p volume
    ${install}

    # install only extracts; make sure every chromosome file matches the registered checksums
    python -c 'from SigProfilerMatrixGenerator.scripts import reference_genome_manager as rgm; assert rgm.ReferenceGenomeManager("volume").is_genome_installed("${genome}"), "${genome} failed the SigProfilerMatrixGenerator checksum verification"'
    """

    stub:
    """
    mkdir -p volume/tsb/${genome}
    touch volume/tsb/${genome}/1.txt
    """
}
