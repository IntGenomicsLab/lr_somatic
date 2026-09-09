process SIGPROFILER_ASSIGNMENT {
    tag "$meta.id"
    label 'process_low'

    // Conda is not supported: the image installs SigProfilerMatrixGenerator from the fork
    // that adds the CHM13-T2T genome (SigProfilerSuite/SigProfilerMatrixGenerator#250) and
    // SigProfilerAssignment from the fork that adds CHM13-T2T COSMIC signatures. Return to
    // the bioconda/biocontainers releases once both are merged and released upstream.
    container "${(workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer') && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/ljwharbers/sigprofiler-sif:1.3.6-chm13-28a9ce8'
        : 'ghcr.io/ljwharbers/sigprofiler:1.3.6-chm13-28a9ce8'}"

    input:
    tuple val(meta), path(sbs96), path(dbs78), path(id83)   // SigProfilerMatrixGenerator matrices; dbs78/id83 may be []
    val(genome)                                             // genome build for the COSMIC reference signatures, e.g. GRCh38 or CHM13-T2T
    val(cosmic_version)                                     // COSMIC version, e.g. 3.6

    output:
    tuple val(meta), path("*/SBS96")                                              , emit: sbs96     , optional: true
    tuple val(meta), path("*/DBS78")                                              , emit: dbs78     , optional: true
    tuple val(meta), path("*/ID83")                                               , emit: id83      , optional: true
    tuple val(meta), path("*/*/Assignment_Solution/Activities/*_Activities.txt")  , emit: activities, optional: true
    tuple val("${task.process}"), val('sigprofilerassignment'), eval("python -c 'import importlib.metadata as m; print(m.version(\"SigProfilerAssignment\"))'"), topic: versions, emit: versions_sigprofilerassignment

    when:
    task.ext.when == null || task.ext.when

    script:
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SIGPROFILER_ASSIGNMENT does not support Conda. Please use Docker / Singularity / Apptainer instead."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // SBS matrices are collapsed to SBS96 for the fit; DBS78 and ID83 must not be
    // (SigProfilerAssignment refuses collapse_to_SBS96 for non-SBS contexts).
    def fits = [[sbs96, 'SBS96', 'True'], [dbs78, 'DBS78', 'False'], [id83, 'ID83', 'False']]
        .findAll { matrix, _ctx, _collapse -> matrix }
        .collect { matrix, ctx, collapse ->
            """
    if [ "\$(matrix_total ${matrix})" -gt 0 ]; then
        SigProfilerAssignment cosmic_fit \\
            ${matrix} \\
            ${prefix}/${ctx} \\
            --input_type matrix \\
            --genome_build ${genome} \\
            --cosmic_version ${cosmic_version} \\
            --collapse_to_SBS96 ${collapse} \\
            --volume spa_volume \\
            --cpu ${task.cpus} \\
            ${args}
    else
        echo "${ctx}: no mutations in ${matrix}, skipping COSMIC fit" >&2
    fi
    """
        }
        .join('\n')
    """
    # total mutation count of a SigProfilerMatrixGenerator matrix (all samples, all classes)
    matrix_total() {
        awk -F'\\t' 'NR > 1 { for (i = 2; i <= NF; i++) s += \$i } END { printf "%d\\n", s }' "\$1"
    }

    # sigProfilerPlotting caches plot templates in its (read-only) package directory unless
    # this variable points elsewhere; SigProfilerAssignment does not forward --volume to it.
    mkdir -p ${prefix} spa_volume
    export SIGPROFILERPLOTTING_VOLUME="\$PWD/spa_volume"
    ${fits}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/SBS96/Assignment_Solution/Activities
    touch ${prefix}/SBS96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt
    """
}
