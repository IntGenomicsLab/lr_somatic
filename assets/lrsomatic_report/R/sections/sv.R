# Structural variants section (reference implementation of the section-module contract), keyed by caller

register_section(list(
  id = "sv",
  title = "Structural variants",

  locate = function(sample_dir, sample_id) {
    d = sample_dir

    find1 = function(pattern) {
      hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) hits[1] else NULL
    }

    severus_vcf      = find1("^severus_somatic\\.vcf\\.gz$")
    severus_gene_tsv = find1("^SV_filtered_with_gene_annotations\\.tsv$")
    # VEP SV VCF is the usual annotation source; the gene-annotated TSV is the fallback
    severus_vep_vcf  = find1("_SV_VEP\\.vcf\\.gz$")

    list(callers = list(
      severus = list(vcf = severus_vcf, gene_tsv = severus_gene_tsv, vep_vcf = severus_vep_vcf)
    ))
  },

  parse = function(inputs, section_data) {
    tabs = list()
    circ = list(nontrans = data.table(), translocations = data.table())
    annotation_path = NULL

    for (nm in names(inputs$callers)) {
      caller_inputs = inputs$callers[[nm]]

      # One parse feeds both the table and the circos tracks
      records = parse_severus_somatic_records(caller_inputs$vcf)

      if (nrow(records) > 0) {
        # A missing VEP SV VCF only empties the symbol columns; panel matching is on coordinates
        t = build_sv_table_from_vep(records, caller_inputs$vep_vcf)
        if (!is.null(caller_inputs$vep_vcf)) annotation_path = caller_inputs$vep_vcf
      } else {
        # No caller VCF: fall back to the gene-annotated TSV
        t = build_sv_table(parse_severus_gene_tsv(caller_inputs$gene_tsv))
        if (nrow(t) > 0) annotation_path = caller_inputs$gene_tsv
      }

      if (!is.null(t) && nrow(t) > 0) {
        t[, caller := nm]
        tabs[[nm]] = t
      }
      # Circos tracks from the same collapsed records; last-write-wins is a no-op with one caller
      circ = severus_circos_tracks(records)
    }

    tbl = if (length(tabs) > 0) rbindlist(tabs, fill = TRUE) else data.table()

    list(
      table            = tbl,
      circos           = circ,
      annotation_path  = annotation_path,
      annotation_found = !is.null(annotation_path)
    )
  }
))
