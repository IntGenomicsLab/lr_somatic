# Phasing section (WhatsHap stats from qc/whatshap_stats/); these are germline statistics

register_section(list(
  id = "whatshap",
  title = "Phasing",

  locate = function(sample_dir, sample_id) {
    d = sample_dir

    find1 = function(pattern) {
      hits = list.files(d, pattern = pattern, recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) hits[1] else NULL
    }

    # qc/whatshap_stats/ is not tumour/normal-scoped, so a plain recursive match is correct
    list(stats_tsv = find1("_whatshap_stats\\.tsv$"))
  },

  parse = function(inputs, section_data) {
    f = inputs$stats_tsv
    if (is.null(f) || !file.exists(f)) return(NULL)

    dt = tryCatch(
      fread(f, sep = "\t", header = TRUE),
      error = function(e) {
        message("Failed to parse WhatsHap stats: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(dt) || nrow(dt) == 0) return(NULL)

    # The header line is "#sample\tchromosome\t..." — fread keeps the leading "#".
    setnames(dt, sub("^#", "", names(dt)))
    if (!"chromosome" %in% names(dt)) {
      message("WhatsHap stats has no 'chromosome' column; skipping section")
      return(NULL)
    }

    # bp_per_block_sum reads as integer64, which DT renders badly; widen to double
    for (col in names(dt)) {
      if (inherits(dt[[col]], "integer64")) dt[, (col) := as.numeric(get(col))]
    }

    list(
      per_chrom = dt[chromosome != "ALL"],
      all       = if (any(dt$chromosome == "ALL")) as.list(dt[chromosome == "ALL"][1]) else NULL,
      vcf       = if ("file_name" %in% names(dt)) dt$file_name[1] else NA_character_
    )
  }
))
