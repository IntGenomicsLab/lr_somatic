process LRSOMATICREPORT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    // Dependencies only (the tool is vendored at assets/lrsomatic_report); when environment.yml changes rebuild both images with `wave --conda-file modules/local/lrsomaticreport/environment.yml --freeze --await [--singularity]`
    container "${workflow.containerEngine == 'singularity'
        ? 'oras://community.wave.seqera.io/library/r-base_quarto_r-base64enc_r-data.table_pruned:dc62d809aa6fd497'
        : 'community.wave.seqera.io/library/r-base_quarto_r-base64enc_r-data.table_pruned:c1049dbaf31bf178'}"

    input:
    // Every path input is optional (`[]` when skipped); tumor/normal QC stage into separate dirs because a matched pair shares meta.id
    tuple val(meta), path(vep_somatic), path(sv_vep), path(severus_vcf), path(somatic_vcf), path(ascat_files), path(qc_tumor_files, stageAs: 'qc_tumor/*'), path(qc_normal_files, stageAs: 'qc_normal/*'), path(wakhan_files, stageAs: 'wakhan/*')
    path(report_src) // lrsomatic_report source tree (bin/, R/, templates/, assets/)
    // User-supplied gene panel TSVs (`[]` for builtins); the matching `--gene-panel gene_panels/<base>` args are built in conf/modules.config
    path(gene_panels, stageAs: 'gene_panels/*')

    output:
    tuple val(meta), path("*_report.html"), emit: report
    // WARN: Manually update to match the vendored release in assets/lrsomatic_report/VENDORED.md
    tuple val("${task.process}"), val('lrsomatic_report'), val('1.3.2'), topic: versions, emit: versions_lrsomaticreport

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sex = meta.sex ?: 'male'

    // Discovery is recursive and matches on base name, so suffix-distinct files can be linked flat
    def flat_inputs = [vep_somatic, sv_vep, severus_vcf, ascat_files].flatten().findAll { f -> f }
    def link_flat = flat_inputs ? """
    for f in ${flat_inputs.collect { f -> "\"${f}\"" }.join(' ')}; do ln -s "\$PWD/\$f" "sample_dir/\$f"; done
    """ : ''

    // The VAF/depth/phasing source is looked up at a literal path
    def link_somatic = somatic_vcf ? """
    mkdir -p sample_dir/variants/phased
    ln -s "\$PWD/${somatic_vcf}" sample_dir/variants/phased/somatic_smallvariants.vcf.gz
    """ : ''

    """
    # Quarto/Deno write under \$HOME and \$TMPDIR, which clusters may mount read-only
    export HOME=\$PWD
    export TMPDIR=\$PWD/tmp TMP=\$PWD/tmp TEMP=\$PWD/tmp
    mkdir -p "\$TMPDIR"

    # The Wave container doesn't source conda's activation hooks (quarto needs QUARTO_SHARE_PATH); -profile conda may already set CONDA_PREFIX
    export CONDA_PREFIX="\${CONDA_PREFIX:-/opt/conda}"
    for f in "\$CONDA_PREFIX"/etc/conda/activate.d/*.sh; do
        [ -f "\$f" ] && source "\$f"
    done

    mkdir -p sample_dir
    ${link_flat}
    ${link_somatic}

    # Link file by file: R's list.files(recursive = TRUE) does not descend into symlinked dirs
    if [ -d qc_tumor ]; then
        mkdir -p sample_dir/qc/tumor
        for f in qc_tumor/*; do ln -s "\$PWD/\$f" "sample_dir/qc/tumor/\$(basename "\$f")"; done
    fi
    if [ -d qc_normal ]; then
        mkdir -p sample_dir/qc/normal
        for f in qc_normal/*; do ln -s "\$PWD/\$f" "sample_dir/qc/normal/\$(basename "\$f")"; done
    fi

    # Wakhan is addressed by fixed path: sample_dir/wakhan must hold solutions_ranks.tsv, the heatmap and solution_<rank>/
    if [ -d wakhan ]; then
        mkdir -p sample_dir/wakhan
        for f in wakhan/*; do ln -s "\$PWD/\$f" "sample_dir/wakhan/\$(basename "\$f")"; done
    fi

    Rscript "${report_src}/bin/render_report.R" \\
        --sample-dir sample_dir \\
        --sample-id "${prefix}" \\
        --sex "${sex}" \\
        --reference auto \\
        --output "${prefix}_report.html" \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_report.html
    """
}
