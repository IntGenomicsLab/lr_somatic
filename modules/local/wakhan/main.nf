process WAKHAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/wakhan:0.4.3--pyhdfd78af_0':
        'biocontainers/wakhan:0.4.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(tumor_input), path(tumor_index), path(normal_input), path(normal_index), path(vcf), path(breakpoints)
    tuple val(meta2), path(reference)
    path(centromere_bed)

    output:
    tuple val(meta), path("*/*_genes_genome.html")                              , emit: genes_genome_html
    tuple val(meta), path("*/*_genes_genome.pdf")                               , emit: genes_genome_pdf
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints.html")            , emit: breakpoints_html
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints.pdf")             , emit: breakpoints_pdf
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints_subclonal.html")  , emit: breakpoints_subclonal_html
    tuple val(meta), path("*/*_genome_copynumbers_breakpoints_subclonal.pdf")   , emit: breakpoints_subclonal_pdf
    tuple val(meta), path("*/*_genome_copynumbers_details.html")                , emit: copynumbers_details_html
    tuple val(meta), path("*/*_genome_copynumbers_details.pdf")                 , emit: copynumbers_details_pdf
    tuple val(meta), path("*/bed_output/*.bed")                                 , emit: bed_files
    tuple val(meta), path("*/variation_plots/*.html")                           , emit: variation_plots
    tuple val(meta), path("*/vcf_output/*_wakhan_cna_*.vcf")                    , emit: vcf_files
    tuple val(meta), path("*_heatmap_ploidy_purity.html")                       , emit: heatmap_html
    tuple val(meta), path("*_heatmap_ploidy_purity.html.pdf")                   , emit: heatmap_pdf
    tuple val(meta), path("coverage_data/*.csv")                                , emit: coverage_csv
    tuple val(meta), path("coverage_data/*.png")                                , emit: coverage_png
    tuple val(meta), path("coverage_plots/*.html")                              , emit: coverage_plots_html
    tuple val(meta), path("coverage_plots/*.pdf")                               , emit: coverage_plots_pdf
    tuple val(meta), path("phasing_output/*.html")                              , emit: phasing_html
    tuple val(meta), path("phasing_output/*.pdf")                               , emit: phasing_pdf
    tuple val(meta), path("phasing_output/*rephased.vcf.gz")                    , emit: rephased_vcf
    tuple val(meta), path("phasing_output/*rephased.vcf.gz.csi")                , emit: rephased_vcf_index
    tuple val(meta), path("snps_loh_plots/*_genome_snps_ratio_loh.html")        , emit: snps_loh_plot,      optional: true
    tuple val(meta), path("solutions_ranks.tsv")                                , emit: solutions_ranks
    // Whole directories, not the plots inside: every solution's plot has the same basename,
    // and LRSOMATICREPORT resolves them by solution_<rank>/ path
    tuple val(meta), path("solution_*", type: 'dir')                            , emit: solution_dirs,   optional: true
    // WARN: Manually update version information as tool does not provide on CLI
    tuple val("${task.process}"), val('wakhan'), val("0.4.3"), topic: versions, emit: versions_wakhan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def phased_vcf = normal_input ? "--normal-phased-vcf $vcf" : "--tumor-phased-vcf $vcf"
    def centromere = centromere_bed ? "--centromere \$PWD/${centromere_bed}" : ""

    """
    wakhan \\
        all \\
        --target-bam ${tumor_input} \\
        --breakpoints ${breakpoints} \\
        --reference ${reference} \\
        --genome-name ${prefix} \\
        --out-dir-plots . \\
        ${phased_vcf} \\
        ${centromere} \\
        ${args} \\
        --threads ${task.cpus}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    """
}
