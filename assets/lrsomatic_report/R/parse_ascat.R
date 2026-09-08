suppressPackageStartupMessages({
  library(data.table)
})

# Parse ASCAT raw segments (segments_raw.txt)
parse_ascat_segments = function(segments_file) {
  if (is.null(segments_file) || !file.exists(segments_file)) return(NULL)
  dt = fread(segments_file, sep = "\t", header = TRUE)

  # Normalise column names to lowercase
  setnames(dt, tolower(names(dt)))

  # Add chr prefix if missing
  dt[, chr := ensure_chr_prefix(as.character(chr))]

  # Column names after tolower(): naraw, nbraw
  dt[, total_cn := pmin(naraw + nbraw, 4)]
  dt[, major_cn := pmin(naraw, 4)]
  dt[, minor_cn := pmin(nbraw, 4)]

  dt
}

# Parse ASCAT purity/ploidy file
parse_ascat_purityploidy = function(pp_file) {
  if (is.null(pp_file) || !file.exists(pp_file)) return(list(purity = NA_real_, ploidy = NA_real_))
  dt = fread(pp_file, sep = "\t", header = TRUE)
  setnames(dt, tolower(names(dt)))
  list(
    purity = round(as.numeric(dt$aberrantcellfraction[1]), 3),
    ploidy = round(as.numeric(dt$ploidy[1]), 3)
  )
}

# Parse Wakhan's ranked purity/ploidy solutions table (wakhan/solutions_ranks.tsv)
parse_wakhan_solutions = function(tsv_file) {
  if (is.null(tsv_file) || !file.exists(tsv_file)) return(NULL)
  dt = fread(tsv_file, sep = "\t", header = TRUE)
  if (nrow(dt) == 0) return(NULL)
  setorder(dt, solution_rank)
  dt
}

# Locate each solution's genome copy-number plot; try solution_<rank>/ first to avoid the aliased duplicate directory
locate_wakhan_cn_plots = function(wakhan_dir, solutions_dt) {
  if (is.null(wakhan_dir) || is.null(solutions_dt) || nrow(solutions_dt) == 0) return(list())
  out = lapply(seq_len(nrow(solutions_dt)), function(i) {
    row = solutions_dt[i]
    sdir = file.path(wakhan_dir, paste0("solution_", row$solution_rank))
    if (!dir.exists(sdir)) sdir = file.path(wakhan_dir, row$repository_name)
    if (!dir.exists(sdir)) return(NULL)
    hits = list.files(sdir, pattern = "genome_copynumbers_breakpoints\\.html$", full.names = TRUE)
    if (length(hits) == 0) return(NULL)
    list(rank = row$solution_rank, purity = row$cell_purity, ploidy = row$ploidy, plot = hits[1])
  })
  Filter(Negate(is.null), out)
}
