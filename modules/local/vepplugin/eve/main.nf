process VEPPLUGIN_EVE {
    tag "${eve_dir}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f'}"

    input:
    path eve_dir

    output:
    path "eve_merged.vcf.gz{,.tbi}", emit: files
    tuple val("${task.process}"), val('tabix'), eval("tabix -h 2>&1 | grep -oP 'Version:\\s*\\K[^\\s]+'"), topic: versions, emit: versions_tabix

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # Merging thousands of per-protein VCFs needs scratch space rather than /tmp
    export TMPDIR=\${TMPDIR:-\$PWD}

    prepare_vep_plugin_data.sh eve ${eve_dir} . ${args}
    """

    stub:
    """
    touch eve_merged.vcf.gz
    touch eve_merged.vcf.gz.tbi
    """
}
