suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# INFO/CALLER values that denote a somatic caller (the merged VEP VCF is mostly germline records); ClairS is "clairs" or "clairs-to"
SOMATIC_CALLERS = c("clairs", "clairs-to", "clairsto", "deepsomatic")

# Derive dbsnp/cosmic columns from VEP's Existing_variation list
derive_dbsnp_cosmic = function(dt) {
  dt[, dbsnp := sub("(rs[0-9]+).*", "\\1", existing)]
  dt[!grepl("^rs", dbsnp, perl = TRUE), dbsnp := NA_character_]

  dt[, cosmic := sub(".*(COS[VM][0-9]+).*", "\\1", existing)]
  dt[!grepl("^COS", cosmic, perl = TRUE), cosmic := NA_character_]
  dt
}

# Dispatch on file contents: VEP text output (#Uploaded_variation header) vs VCF with CSQ (#CHROM header)
parse_vep = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  con = gzfile(vep_file, "rb")
  is_vcf = FALSE
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
    if (startsWith(line, "#CHROM")) { is_vcf = TRUE; break }
  }
  close(con)

  if (is_vcf) parse_vep_vcf(vep_file) else parse_vep_text(vep_file)
}

# Parse VEP default text output (tab-delimited, not a VCF); one row per consequence per variant
parse_vep_text = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Count meta-lines (start with ##) to find the column-header line
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#Uploaded_variation")) break
    skip_n = skip_n + 1L
  }
  close(con)

  dt = tryCatch(
    fread(vep_file, skip = skip_n, sep = "\t", header = TRUE,
          col.names = function(x) gsub("^#", "", x)),
    error = function(e) {
      message("Failed to parse VEP file: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  setnames(dt, old = "Uploaded_variation", new = "variant_id", skip_absent = TRUE)
  setnames(dt, old = "Gene",              new = "gene_id",     skip_absent = TRUE)
  setnames(dt, old = "Consequence",       new = "consequence", skip_absent = TRUE)

  # Coordinates and alleles both from variant_id where it has the canonical shape; Location is off by one for dash-form insertions, so it is only the fallback
  vid = "^.+_[0-9]+_[^_]+/[^_]+$"
  dt[, from_vid := grepl(vid, variant_id)]

  dt[, chrom := ifelse(from_vid, sub("_[0-9]+_[^_]+$", "", variant_id),
                                 sub(":.*", "", Location))]
  dt[, pos   := as.integer(ifelse(from_vid, sub(".*_([0-9]+)_[^_]+$", "\\1", variant_id),
                                            sub(".*:(\\d+).*", "\\1", Location)))]
  dt[, chrom := ensure_chr_prefix(chrom)]

  dt[, ref := sub(".*_([^/]+)/.*", "\\1", variant_id)]
  dt[, alt := sub(".*/", "", variant_id)]

  # Parse VEP Extra key=value field
  dt[, symbol   := extract_extra_key(Extra, "SYMBOL")]
  dt[, impact   := extract_extra_key(Extra, "IMPACT")]
  dt[, existing := extract_extra_key(Extra, "Existing_variation")]
  dt[, sift     := extract_extra_key(Extra, "SIFT")]
  dt[, polyphen := extract_extra_key(Extra, "PolyPhen")]
  dt[, hgvsp    := extract_extra_key(Extra, "HGVSp")]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt = derive_dbsnp_cosmic(dt)

  # No per-variant caller in the text format; kept for contract parity with parse_vep_vcf()
  dt[, caller := NA_character_]

  # Record identity for parity with parse_vep_vcf(): VEP's own Uploaded_variation name, never a caller ID
  dt[, id := variant_id]

  # chrom/pos/ref/alt are in VEP notation (indels shifted, possibly dash-form); see coord_space in build_variant_table()
  dt[, coord_space := "vep"]

  dt
}

# Parse a VCF with VEP CSQ annotation; same column contract as parse_vep_text()
parse_vep_vcf = function(vep_file) {
  if (is.null(vep_file) || !file.exists(vep_file)) return(NULL)

  # Skip header to #CHROM, capturing the CSQ field order and whether a CALLER tag exists
  con = gzfile(vep_file, "rb")
  skip_n = 0L
  csq_format = NULL
  has_caller_info = FALSE
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "##INFO=<ID=CSQ")) {
      m = regmatches(line, regexpr("Format: [^\"]+", line))
      if (length(m) > 0) csq_format = strsplit(sub("^Format: ", "", m), "|", fixed = TRUE)[[1]]
    }
    if (startsWith(line, "##INFO=<ID=CALLER")) has_caller_info = TRUE
    if (startsWith(line, "#CHROM")) break
    skip_n = skip_n + 1L
  }
  close(con)

  if (is.null(csq_format)) {
    message("Failed to parse VEP VCF: no CSQ Format found in header")
    return(NULL)
  }

  dt = tryCatch(
    fread(vep_file, skip = skip_n + 1L, sep = "\t", header = FALSE, select = 1:8,
          col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")),
    error = function(e) {
      message("Failed to parse VEP VCF: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  # Drop records the callers themselves rejected (RefCall / LowQual / GERMLINE).
  dt = dt[FILTER %in% c("PASS", ".")]
  if (nrow(dt) == 0) return(NULL)

  # Fixed-string splits: much faster than sub() over millions of rows
  dt[, CSQ := tstrsplit(INFO, "CSQ=", fixed = TRUE, keep = 2L)[[1]]]
  dt[, CSQ := tstrsplit(CSQ,  ";",    fixed = TRUE, keep = 1L)[[1]]]

  # Keep somatic-caller records only; skipped for older VCFs with no CALLER tag at all
  if (has_caller_info) {
    dt[, caller := tstrsplit(INFO, "CALLER=", fixed = TRUE, keep = 2L)[[1]]]
    dt[, caller := tstrsplit(caller, ";", fixed = TRUE, keep = 1L)[[1]]]
    dt = dt[caller %in% SOMATIC_CALLERS]
    if (nrow(dt) == 0) return(NULL)
  } else {
    dt[, caller := NA_character_]
  }

  # Collapse multi-caller duplicates to one row per variant before CSQ expansion; ID rides along for the SV join
  dt = dt[, .(CSQ = CSQ[1], id = ID[1],
              caller = paste(sort(unique(caller)), collapse = ",")),
          by = .(CHROM, POS, REF, ALT)]
  dt[id == ".", id := NA_character_]

  # One row per gene/transcript annotation (comma-separated CSQ entries)
  dt_long = dt[, .(csq_entry = unlist(strsplit(CSQ, ",", fixed = TRUE))),
               by = .(CHROM, POS, REF, ALT, id, caller)]
  if (nrow(dt_long) == 0) return(NULL)

  # Split each entry on "|", keeping only the fields used below
  need = c("Consequence", "IMPACT", "SYMBOL", "Gene", "HGVSp",
           "Existing_variation", "SIFT", "PolyPhen")
  keep_idx = match(need, csq_format)
  ok = !is.na(keep_idx)
  parts = tstrsplit(dt_long$csq_entry, "|", fixed = TRUE, fill = NA_character_,
                    keep = keep_idx[ok])
  names(parts) = need[ok]

  get_field = function(name) {
    if (is.null(parts[[name]])) rep(NA_character_, nrow(dt_long)) else parts[[name]]
  }

  dt_long[, consequence := gsub("&", ",", get_field("Consequence"))]
  dt_long[, impact      := get_field("IMPACT")]
  dt_long[, symbol      := get_field("SYMBOL")]
  dt_long[, gene_id     := get_field("Gene")]
  dt_long[, hgvsp       := get_field("HGVSp")]
  dt_long[, existing    := get_field("Existing_variation")]
  dt_long[, sift        := get_field("SIFT")]
  dt_long[, polyphen    := get_field("PolyPhen")]

  # Blank CSQ fields are "", not NA; normalise for parity with parse_vep_text()
  for (col in c("symbol", "impact", "hgvsp", "existing", "gene_id", "sift", "polyphen",
                "caller")) {
    dt_long[get(col) == "", (col) := NA_character_]
  }

  dt_long[, chrom := ensure_chr_prefix(CHROM)]
  dt_long[, pos   := as.integer(POS)]
  dt_long[, ref   := REF]
  dt_long[, alt   := ALT]

  # These are the VCF's own POS/REF/ALT (not VEP-shifted); see build_variant_table()
  dt_long[, coord_space := "vcf"]

  # dbSNP / COSMIC IDs, derived from Existing_variation
  dt_long = derive_dbsnp_cosmic(dt_long)

  dt_long[, .(chrom, pos, ref, alt, id, symbol, gene_id, consequence, impact, hgvsp,
              existing, dbsnp, cosmic, sift, polyphen, caller, coord_space)]
}

# Pick the sample column: one is unambiguous; else the report's sample, then the first non-normal; the first column is a flagged guess
pick_sample_column = function(samples, sample_id = NULL, vcf_file = "") {
  if (length(samples) == 1L) return(1L)

  if (!is.null(sample_id) && nzchar(sample_id)) {
    i = match(sample_id, samples)
    if (!is.na(i)) return(i)
  }
  not_normal = which(!grepl("^(normal|blood|germline|ref)", samples, ignore.case = TRUE))
  if (length(not_normal) > 0) {
    if (length(not_normal) > 1L || is.null(sample_id))
      warning("VCF '", vcf_file, "' has samples [", paste(samples, collapse = ", "),
              "]; reading '", samples[not_normal[1]], "'.")
    return(not_normal[1])
  }
  warning("VCF '", vcf_file, "' has samples [", paste(samples, collapse = ", "),
          "] and none looks like a tumour; reading '", samples[1], "'.")
  1L
}

# Allele fraction from whichever FORMAT tag is present: LRSomatic renames AF/VAF per caller (STANDARDIZE_AF); AD is the last resort
allele_fraction = function(field, fields) {
  for (tag in c("AF", "VAF")) {
    if (tag %in% fields) {
      # Multi-allelic: the first value pairs with the first ALT, which is all this table joins on
      v = suppressWarnings(as.numeric(sub(",.*$", "", field(tag))))
      if (any(!is.na(v))) return(v)
    }
  }
  if ("AD" %in% fields) {
    ad = strsplit(field("AD"), ",", fixed = TRUE)
    return(vapply(ad, function(x) {
      n = suppressWarnings(as.numeric(x))
      if (length(n) < 2L || anyNA(n[1:2]) || sum(n[1:2]) == 0) return(NA_real_)
      n[2] / sum(n[1:2])
    }, numeric(1)))
  }
  rep(NA_real_, length(field("GT")))
}

# Parse raw caller VCF(s) into chrom, pos, ref, alt, vaf, dp, gt, ps, caller; several paths are stacked (ClairS splits snvs/indels). `sample_id` selects the column by name so a matched VCF never reports the normal's values
parse_caller_vcf = function(vcf_file, caller_name = "unknown", sample_id = NULL) {
  if (is.null(vcf_file)) return(NULL)
  if (length(vcf_file) > 1) {
    parts = lapply(vcf_file, parse_caller_vcf, caller_name = caller_name,
                   sample_id = sample_id)
    parts = parts[!vapply(parts, is.null, logical(1))]
    return(if (length(parts) > 0) rbindlist(parts) else NULL)
  }
  if (!file.exists(vcf_file)) return(NULL)

  # Count header lines, keeping the #CHROM line itself — it carries the sample names.
  con = gzfile(vcf_file, "rb")
  skip_n = 0L
  chrom_line = NULL
  repeat {
    line = readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "#CHROM")) { chrom_line = line; break }
    skip_n = skip_n + 1L
  }
  close(con)
  if (is.null(chrom_line)) return(NULL)

  header  = strsplit(chrom_line, "\t", fixed = TRUE)[[1]]
  samples = if (length(header) > 9L) header[10:length(header)] else character(0)
  # Sites-only VCF: no genotypes, return NULL
  if (length(samples) == 0) return(NULL)

  sample_col = pick_sample_column(samples, sample_id, vcf_file)

  # 1:9 are the fixed VCF columns; the chosen sample sits at 9 + its index.
  col_names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                "SAMPLE1")
  dt = tryCatch(
    fread(vcf_file, skip = skip_n + 1L, sep = "\t", header = FALSE,
          select = c(1:9, 9L + sample_col), col.names = col_names),
    error = function(e) {
      warning("Could not read VCF '", vcf_file, "': ", conditionMessage(e))
      NULL
    })
  if (is.null(dt) || nrow(dt) == 0) return(NULL)

  dt[, CHROM := ensure_chr_prefix(CHROM)]

  # Extract AF/DP/GT/PS from FORMAT + SAMPLE, format-group by format-group
  fmt_groups = unique(dt$FORMAT)
  vaf_list = rep(NA_real_,      nrow(dt))
  dp_list  = rep(NA_integer_,   nrow(dt))
  gt_list  = rep(NA_character_, nrow(dt))
  ps_list  = rep(NA_character_, nrow(dt))

  for (fmt in fmt_groups) {
    idx_rows = which(dt$FORMAT == fmt)
    fields  = strsplit(fmt, ":", fixed = TRUE)[[1]]
    split_s = strsplit(dt$SAMPLE1[idx_rows], ":", fixed = TRUE)

    # One FORMAT field, by name, across this group's rows
    field = function(name) {
      i = match(name, fields)
      if (is.na(i)) return(rep(NA_character_, length(split_s)))
      vapply(split_s, function(x) if (length(x) >= i) x[i] else NA_character_,
             character(1))
    }

    vaf_list[idx_rows] = allele_fraction(field, fields)
    dp_list[idx_rows]  = suppressWarnings(as.integer(field("DP")))
    gt_list[idx_rows]  = field("GT")
    ps_list[idx_rows]  = field("PS")
  }

  # Blank the unphased placeholders ("." PS, "/" GT) so cells render empty
  ps_list[!is.na(ps_list) & ps_list == "."] = NA_character_
  gt_list[!is.na(gt_list) & gt_list %in% c(".", "./.")] = NA_character_

  data.table(chrom = dt$CHROM, pos = dt$POS, ref = dt$REF, alt = dt$ALT,
             vaf = vaf_list, dp = dp_list, gt = gt_list, ps = ps_list,
             caller = caller_name)
}

# Header-only provenance for the VCF(s) supplying VAF/DP/GT/PS (after a consensus merge the FORMAT fields need not match the `callers` column); NULL when there is nothing to describe
vaf_provenance = function(vcf_files, sample_dir = NULL, sample_id = NULL) {
  if (is.null(vcf_files) || length(vcf_files) == 0) return(NULL)
  vcf_files = vcf_files[!is.na(vcf_files) & nzchar(vcf_files)]
  if (length(vcf_files) == 0) return(NULL)

  # Paths read better relative to the sample directory the caller passed in.
  rel = function(p) {
    if (is.null(sample_dir) || !nzchar(sample_dir)) return(p)
    root = sub("/+$", "", normalizePath(sample_dir, mustWork = FALSE))
    full = normalizePath(p, mustWork = FALSE)
    if (startsWith(full, paste0(root, "/"))) substring(full, nchar(root) + 2L) else p
  }

  # "##source=..." is optional colour; many VCFs carry none
  read_sources = function(p) {
    if (!file.exists(p)) return(character(0))
    con = gzfile(p, "rb")
    on.exit(close(con), add = TRUE)
    out = character(0)
    repeat {
      line = readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) break
      if (!startsWith(line, "##")) break   # #CHROM or a malformed header ends the scan
      if (startsWith(line, "##source=")) out = c(out, sub("^##source=", "", line))
    }
    out
  }

  # Which sample column parse_caller_vcf() read; named in the footnote because a wrong pick is otherwise invisible
  read_sample = function(p) {
    if (!file.exists(p)) return(NA_character_)
    con = gzfile(p, "rb")
    on.exit(close(con), add = TRUE)
    repeat {
      line = readLines(con, n = 1, warn = FALSE)
      if (length(line) == 0) return(NA_character_)
      if (startsWith(line, "#CHROM")) {
        header = strsplit(line, "\t", fixed = TRUE)[[1]]
        if (length(header) <= 9L) return(NA_character_)
        samples = header[10:length(header)]
        # Only worth reporting when there was a choice to get wrong.
        if (length(samples) == 1L) return(NA_character_)
        return(samples[suppressWarnings(pick_sample_column(samples, sample_id, p))])
      }
      if (!startsWith(line, "##")) return(NA_character_)
    }
  }
  chosen = unique(unlist(lapply(vcf_files, read_sample)))
  chosen = chosen[!is.na(chosen)]

  list(
    paths   = unname(vapply(vcf_files, rel, character(1))),
    sources = unique(unlist(lapply(vcf_files, read_sources))),
    sample  = chosen
  )
}

# Canonical variant key joining VEP rows to VCF records: VEP reports indels at anchor + 1 with raw or dash-trimmed alleles; all forms reconcile as trimmed alleles at anchor + 1, SNVs/MNVs verbatim. `space` is "vcf" or "vep"
variant_key = function(chrom, pos, ref, alt, space = c("vcf", "vep")) {
  space = match.arg(space)
  ref = toupper(as.character(ref)); alt = toupper(as.character(alt))
  pos = as.integer(pos)

  trim = function(x) { t = substr(x, 2L, nchar(x)); ifelse(t == "", "-", t) }

  dash     = ref == "-" | alt == "-"   # already trimmed by VEP
  is_indel = dash | nchar(ref) != nchar(alt)

  # Raw allele pairs still need the anchor base dropped; dash forms are already trimmed.
  need_trim = is_indel & !dash

  # Only the VCF side needs shifting — VEP has already done it.
  key_pos = ifelse(space == "vcf" & is_indel, pos + 1L, pos)
  key_ref = ifelse(need_trim, trim(ref), ref)
  key_alt = ifelse(need_trim, trim(alt), alt)

  paste(chrom, key_pos, key_ref, key_alt, sep = "|")
}

# Classify SNV into 6 SBS mutation categories (C/T-ref normalised)
classify_mut = function(ref, alt) {
  comp = c(A = "T", T = "A", C = "G", G = "C")
  ref = toupper(ref); alt = toupper(alt)
  use_comp = !(ref %in% c("C", "T"))
  norm_ref = ifelse(use_comp, comp[ref], ref)
  norm_alt = ifelse(use_comp, comp[alt], alt)
  paste0(norm_ref, ">", norm_alt)
}

# Small-variant display table: VEP rows with VAF/depth/phasing joined from the annotated VCF; gene_panel filters symbols (NULL = all)
build_variant_table = function(vep_data, vaf_data, gene_panel = NULL) {
  if (is.null(vep_data) || nrow(vep_data) == 0) return(NULL)

  # Impact ranking for deduplication
  impact_rank = c(HIGH = 1L, MODERATE = 2L, LOW = 3L, MODIFIER = 4L)
  vep_data[, impact_rank := impact_rank[impact]]
  vep_data[is.na(impact_rank), impact_rank := 5L]

  # Filter to gene panel (by gene symbol or Ensembl ID fallback)
  if (!is.null(gene_panel)) {
    if (length(gene_panel) > 0) {
      vep_data = vep_data[symbol %in% gene_panel | gene_id %in% gene_panel]
    } else {
      vep_data = vep_data[FALSE]  # Empty panel → empty result
    }
  }
  if (nrow(vep_data) == 0) return(data.table())

  # Keep best consequence per variant×gene (lowest impact rank)
  key_cols = c("chrom", "pos", "ref", "alt", "symbol")
  setorder(vep_data, impact_rank)
  vep_data = unique(vep_data, by = key_cols)

  # Join on the canonical key; the coordinate space is declared by each parser (text = vep, CSQ = vcf), never assumed
  vep_space = if ("coord_space" %in% names(vep_data)) unique(vep_data$coord_space) else "vep"
  stopifnot(length(vep_space) == 1L, vep_space %in% c("vep", "vcf"))
  vep_data[, join_key := variant_key(chrom, pos, ref, alt, space = vep_space)]

  if (!is.null(vaf_data) && nrow(vaf_data) > 0) {
    vdt = vaf_data[, .(join_key = variant_key(chrom, pos, ref, alt, space = "vcf"),
                       vaf, dp, gt, ps)]
    vdt = unique(vdt, by = "join_key")
    vep_data = merge(vep_data, vdt, by = "join_key", all.x = TRUE)
  } else {
    vep_data[, `:=`(vaf = NA_real_, dp = NA_integer_,
                    gt  = NA_character_, ps = NA_character_)]
  }

  # Which caller reported each variant (INFO/CALLER); empty on the VEP text path
  if ("caller" %in% names(vep_data) && any(!is.na(vep_data$caller))) {
    vep_data[, callers := caller]
  } else {
    vep_data[, callers := ""]
  }

  # Mutation category for SNVs
  vep_data[nchar(ref) == 1 & nchar(alt) == 1,
           mut_cat := classify_mut(ref, alt)]

  display_cols = c("symbol", "chrom", "pos", "ref", "alt",
                   "consequence", "impact", "hgvsp",
                   "vaf", "dp", "gt", "ps",
                   "callers", "cosmic", "dbsnp", "sift", "polyphen")
  display_cols = display_cols[display_cols %in% names(vep_data)]
  vep_data[, ..display_cols]
}
