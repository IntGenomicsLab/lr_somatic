//
// Subworkflow with functionality specific to the IntGenomicsLab/lrsomatic pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    _monochrome_logs  // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    _input            //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command,
        false
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    // Parse the input samplesheet CSV and build a per-sample BAM channel
    // Each samplesheet row describes one tumor (+ optional normal) sample
    // Columns: sample_id, bam_tumor, bam_normal, method, sex, fiber,
    //          clair3_model, clairSTO_model, clairS_model, tumor_replicate, normal_replicate
    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        // Step 1: build a combined meta map from the samplesheet columns
        // paired_data = true if a normal BAM is present; false for tumor-only
        .map { meta, bam_tumor, bam_normal, method, sex, fiber, clair3_model, clairSTO_model, clairS_model, tumor_replicate, normal_replicate ->
            def real_clair3_model = (clair3_model == null ) ? null : clair3_model
            def real_clairS_model = (clairS_model == null ) ? null : clairS_model
            def real_clairSTO_model = (clairSTO_model == null ) ? null : clairSTO_model
            def paired_data = bam_normal ? true : false
            def meta_info = meta + [ paired_data: paired_data,
                                     platform: method,         // 'ont' or 'pb'
                                     sex: sex,                 // 'XX', 'XY', or null (for ASCAT)
                                     fiber: fiber,             // 'y' or 'n' (fiber-seq data flag)
                                     clair3_model: real_clair3_model,
                                     clairS_model: real_clairS_model,
                                     clairSTO_model: real_clairSTO_model,
                                     tumor_replicate: tumor_replicate,
                                     normal_replicate: normal_replicate]
            return [ meta_info, [ bam_tumor ], [ bam_normal ?: [] ] ]
        }
        // Flatten BAM lists (handles multi-run entries where bam_tumor/bam_normal are lists)
        .map { meta, bam_tumor, bam_normal ->
           [ meta, bam_tumor.flatten(), bam_normal.flatten() ]
        }
        // Step 2: split each row into separate tumor and normal items
        // flatMap emits 1 item (tumor-only) or 2 items (tumor + normal) per samplesheet row
        // Each item gets type='tumor' or type='normal' and the appropriate replicate ID
        .flatMap { meta, tumor_bam, normal_bam ->
            def meta_tumor = meta.clone()
            meta_tumor.type = 'tumor'
            meta_tumor.replicate = meta_tumor.tumor_replicate
            meta_tumor = meta_tumor.subMap('id',
                                           'paired_data',
                                           'type',
                                           'platform',
                                           'sex',
                                           'fiber',
                                           'clair3_model',
                                           'clairS_model',
                                           'clairSTO_model',
                                           'replicate')
            def result = [[meta_tumor, tumor_bam]]
            // result so far: [[meta_tumor, [tumor_bam_path...]]]

            if (normal_bam) {
                def meta_normal = meta.clone()
                meta_normal.type = 'normal'
                meta_normal.replicate = meta_normal.normal_replicate
                meta_normal = meta_normal.subMap('id',
                                                 'paired_data',
                                                 'type',
                                                 'platform',
                                                 'sex',
                                                 'fiber',
                                                 'clair3_model',
                                                 'clairS_model',
                                                 'clairSTO_model',
                                                 'replicate')
                result << [meta_normal, normal_bam]
                // result now: [[meta_tumor, [tumor_bams]], [meta_normal, [normal_bams]]]
            }

            return result
        }
        .set { ch_samplesheet }

    // Count replicates per sample+type and embed the count in meta as n_replicates.
    // This allows downstream groupTuple() to use groupKey() for eager per-sample release
    // instead of waiting for ALL samples to finish (global synchronization barrier).
    // This groupTuple is safe: the source is a fully-materialized list from
    // samplesheetToList(), so the channel closes immediately without blocking any process.
    ch_samplesheet
        .map { meta, bams -> [[meta.id, meta.type], meta, bams] }
        .groupTuple(by: 0)
        .flatMap { key, metas, bams_list ->
            def n = metas.size()
            [metas, bams_list].transpose().collect { m, b ->
                [m + [n_replicates: n], b]
            }
        }
        .set { ch_samplesheet }

    // ch_samplesheet: [meta, [bam...]]
    //   meta fields: id, paired_data, type ('tumor'|'normal'), platform ('ont'|'pb'),
    //                sex, fiber ('y'|'n'), clair3_model, clairS_model, clairSTO_model,
    //                replicate, n_replicates
    //   paired_data: true for both items in a T/N pair (same value for tumor AND normal rows)
    //   n_replicates: total number of replicates for this sample+type combination
    //   bam: list of paths (multiple runs for same sample remain as a list until SAMTOOLS_CAT)
    //
    // NOTE: tumor-only rows emit ONE item (type='tumor', paired_data=false)
    //       paired rows emit TWO items — tumor (paired_data=true) + normal (paired_data=true)
    //       Both share the same 'id' to allow downstream joins

    emit:
    samplesheet = ch_samplesheet  // [meta, [bam...]]  -- see channel structure above
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    validateReportGenePanels()
}

//
// Split --report_gene_panel into panel tokens (mirrored in conf/modules.config)
//
def reportGenePanelTokens(panel_spec) {
    if (!panel_spec) {
        return []
    }
    return panel_spec.toString().split(',').collect { it.trim() }.findAll { it }
}

//
// Does a --report_gene_panel entry name a file rather than a builtin? Textual because conf/modules.config makes the same call without file()
//
def reportGenePanelIsFile(tok) {
    return tok.contains('/') || tok.toLowerCase().endsWith('.tsv')
}

//
// Builtin panel names bundled with the report tool, reference suffix dropped
//
def reportBuiltinGenePanels() {
    def gene_lists_dir = file("${params.report_src}/assets/gene_lists")
    if (!gene_lists_dir.exists()) {
        return []
    }
    return gene_lists_dir
        .list()
        .findAll { it.endsWith('.tsv') }
        .collect { it.replaceFirst(/(\.(hg38|t2t))?\.tsv$/, '') }
        .unique()
        .sort()
}

//
// Validate --report_gene_panel at launch so a typo fails before alignment and calling run
//
def validateReportGenePanels() {
    if (params.skip_report) {
        return
    }
    def tokens = reportGenePanelTokens(params.report_gene_panel)
    if (!tokens) {
        return
    }

    // "none" means unfiltered, so combining it with a real panel is contradictory
    if (tokens.size() > 1 && tokens.any { it.toLowerCase() == 'none' }) {
        error("--report_gene_panel: 'none' means unfiltered and cannot be combined with other panels, got '${params.report_gene_panel}'. Drop the 'none'.")
    }

    def named = tokens.findAll { tok -> tok.toLowerCase() != 'none' && !reportGenePanelIsFile(tok) }
    def panel_files = tokens.findAll { tok -> reportGenePanelIsFile(tok) }

    def missing = panel_files.findAll { tok -> !file(tok).exists() }
    if (missing) {
        error("--report_gene_panel: panel file not found: '${missing.join("', '")}'.")
    }

    def builtins = reportBuiltinGenePanels()
    def unknown = named.findAll { tok -> !builtins.contains(tok) }
    if (unknown) {
        error("--report_gene_panel: '${unknown.join("', '")}' is not a builtin panel. Builtin panels: ${builtins ? builtins.join(', ') : '<none found>'}. To use a panel file give its path, or a name ending in '.tsv'; use 'none' for no filtering.")
    }

    // Panel files are staged side by side into gene_panels/, so equal base names collide
    def duplicates = panel_files
        .collect { tok -> file(tok).name }
        .countBy { name -> name }
        .findAll { _name, count -> count > 1 }
        .keySet()
    if (duplicates) {
        error("--report_gene_panel: panel files sharing a base name cannot be used together ('${duplicates.join("', '")}'). Rename one of them.")
    }
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, bams) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], bams ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "MultiQC (Ewels et al. 2016),",
            "Samtools (Li et al. 2009),",
            "Mosdepth (Pedersen and Quinlan 2018)."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>",
            "<li>Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., ... & Durbin, R. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25(16), 2078-2079. doi: 10.1093/bioinformatics/btp352</li>",
            "<li>Pedersen, B. S., & Quinlan, A. R. (2018). Mosdepth: quick coverage calculation for genomes and exomes. Bioinformatics, 34(5), 867-868. doi: 10.1093/bioinformatics/btx699</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
