suppressPackageStartupMessages({
  library(data.table)
})

# Parse mosdepth summary; `keep_chroms` avoids a bare "^chr" match (admits decoys, empty on chr-less references)
parse_mosdepth_summary = function(summary_file, keep_chroms = NULL) {
  if (is.null(summary_file) || !file.exists(summary_file)) {
    return(list(mean_depth = NA_real_, table = data.table()))
  }
  dt = fread(summary_file, sep = "\t", header = TRUE)
  setnames(dt, tolower(names(dt)))
  total_row = dt[chrom == "total"]
  mean_depth = if (nrow(total_row) > 0) total_row$mean[1] else NA_real_

  # Keep per-chromosome rows (exclude region-level and total)
  chr_rows = dt[chrom != "total" & !grepl("_region", chrom)]
  if (!is.null(keep_chroms) && length(keep_chroms) > 0) {
    chr_rows = chr_rows[ensure_chr_prefix(chrom) %in% ensure_chr_prefix(keep_chroms)]
  } else {
    chr_rows = chr_rows[grepl("^chr", chrom)]
  }
  total_length = if (nrow(total_row) > 0) total_row$length[1] else NA_real_
  total_bases  = if (nrow(total_row) > 0) total_row$bases[1]  else NA_real_
  list(mean_depth = round(mean_depth, 2), total_length = total_length, total_bases = total_bases, table = chr_rows)
}

# Parse mosdepth global distribution (*.mosdepth.global.dist.txt)
parse_mosdepth_dist = function(dist_file) {
  if (is.null(dist_file) || !file.exists(dist_file)) return(NULL)
  dt = fread(dist_file, sep = "\t", header = FALSE,
             col.names = c("chrom", "coverage", "fraction"))
  dt
}

# Parse cramino alignment report
parse_cramino = function(cramino_file) {
  if (is.null(cramino_file) || !file.exists(cramino_file)) {
    return(list(n50 = NA_real_, yield_gb = NA_real_,
                mapped_pct = NA_real_, n_reads = NA_integer_))
  }
  lines = readLines(cramino_file, warn = FALSE)
  get_val = function(pattern) {
    hit = grep(pattern, lines, value = TRUE, ignore.case = TRUE)
    if (length(hit) == 0) return(NA_character_)
    trimws(sub(paste0(".*", pattern, "\\s*"), "", hit[1], ignore.case = TRUE))
  }

  # Cramino outputs key\tvalue pairs
  dt = tryCatch(
    fread(cramino_file, sep = "\t", header = FALSE, col.names = c("key", "value"), fill = TRUE),
    error = function(e) NULL
  )
  if (is.null(dt)) return(list(n50 = NA_real_, yield_gb = NA_real_,
                               mapped_pct = NA_real_, n_reads = NA_integer_))

  get_field = function(pattern) {
    row = dt[grepl(pattern, key, ignore.case = TRUE)]
    if (nrow(row) == 0) NA_character_ else as.character(row$value[1])
  }

  list(
    n50        = suppressWarnings(as.numeric(get_field("N50"))),
    yield_gb   = suppressWarnings(as.numeric(get_field("Yield"))),
    mapped_pct = suppressWarnings(as.numeric(sub("%", "", get_field("% from total")))),
    n_reads    = suppressWarnings(as.integer(get_field("Number of reads")))
  )
}

# Parse samtools flagstat
parse_flagstat = function(flagstat_file) {
  if (is.null(flagstat_file) || !file.exists(flagstat_file)) return(list())
  lines = readLines(flagstat_file, warn = FALSE)
  out = list()
  for (line in lines) {
    count = suppressWarnings(as.integer(sub(" .*", "", trimws(line))))
    if (grepl("in total", line))    out$total    = count
    if (grepl("mapped \\(", line))  out$mapped   = count
    if (grepl("paired in seq", line)) out$paired = count
    if (grepl("secondary", line))   out$secondary = count
    if (grepl("supplementary", line)) out$supplementary = count
    if (grepl("duplicate", line))   out$duplicate = count
  }
  out
}

# Parse samtools stats (*.stats), SN summary lines only
parse_samtools_stats = function(stats_file) {
  if (is.null(stats_file) || !file.exists(stats_file)) return(NULL)
  lines = readLines(stats_file, warn = FALSE)
  sn = lines[startsWith(lines, "SN\t")]
  get_sn = function(key) {
    hit = grep(paste0("^SN\t", key, ":\t"), sn, value = TRUE)
    if (length(hit) == 0) return(NA_real_)
    suppressWarnings(as.numeric(trimws(sub(paste0("^SN\t", key, ":\t([^\t#]+).*"), "\\1", hit[1]))))
  }
  reads_total  = get_sn("raw total sequences")
  reads_mapped = get_sn("reads mapped")
  mapped_pct   = if (!is.na(reads_total) && reads_total > 0)
    round(reads_mapped / reads_total * 100, 2) else NA_real_
  list(
    reads_total  = reads_total,
    reads_mapped = reads_mapped,
    mapped_pct   = mapped_pct,
    total_length = get_sn("total length"),
    bases_mapped = get_sn("bases mapped \\(cigar\\)"),
    error_rate   = get_sn("error rate"),
    avg_length   = get_sn("average length"),
    max_length   = get_sn("maximum length"),
    avg_quality  = get_sn("average quality")
  )
}
