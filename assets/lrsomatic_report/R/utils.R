suppressPackageStartupMessages({
  library(data.table)
})

# Type-stable: fread() may read CHROM as integer and ifelse() on all-NA returns logical
ensure_chr_prefix = function(x) {
  x   = as.character(x)
  out = as.character(ifelse(startsWith(x, "chr"), x, paste0("chr", x)))
  out[is.na(x)] = NA_character_
  out
}

strip_chr_prefix = function(x) {
  sub("^chr", "", x)
}

# Parse VEP "Extra" key=value semicolon-delimited field into a named character vector
parse_extra_kv = function(extra_string) {
  if (is.na(extra_string) || extra_string == "" || extra_string == "-") return(character(0))
  pairs = strsplit(extra_string, ";", fixed = TRUE)[[1]]
  kv = strsplit(pairs, "=", fixed = TRUE)
  keys = vapply(kv, `[`, character(1), 1)
  vals = vapply(kv, function(x) if (length(x) >= 2) paste(x[-1], collapse = "=") else "", character(1))
  setNames(vals, keys)
}

# Vectorised: extract one key from VEP Extra column for each row
extract_extra_key = function(extra_vec, key) {
  vapply(extra_vec, function(x) {
    kv = parse_extra_kv(x)
    if (key %in% names(kv)) kv[[key]] else NA_character_
  }, character(1), USE.NAMES = FALSE)
}

# ---- Gene panels: plain lists (they round-trip through Quarto execute_params) of name, path, reference, has_coords, genes and, when has_coords, parallel chrom/start/end/interval_gene ----

# Canonical reference names; a coordinate panel declares one and it is checked against the render
normalise_reference_name = function(x) {
  if (is.null(x) || length(x) != 1 || is.na(x) || !nzchar(trimws(x))) return(NA_character_)
  x = tolower(trimws(x))
  if (x %in% c("t2t", "chm13", "chm13v2", "chm13v2.0", "t2t-chm13")) return("t2t")
  if (x %in% c("hg38", "grch38", "hg38-noalt")) return("hg38")
  x
}

# Reference suffix of a builtin panel filename; list(name, reference), reference NA if unsuffixed
.split_panel_filename = function(path) {
  stem  = tools::file_path_sans_ext(basename(path))
  parts = strsplit(stem, ".", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    ref = normalise_reference_name(parts[length(parts)])
    if (!is.na(ref) && ref %in% c("t2t", "hg38"))
      return(list(name = paste(parts[-length(parts)], collapse = "."), reference = ref))
  }
  list(name = stem, reference = NA_character_)
}

# Read the leading "#" comment block of a panel TSV (may declare "# reference: hg38")
.panel_header = function(path) {
  lines = readLines(path, warn = FALSE)
  n_comment = 0L
  declared  = NA_character_
  for (ln in lines) {
    if (!startsWith(ln, "#")) break
    n_comment = n_comment + 1L
    m = regmatches(ln, regexpr("^#\\s*reference\\s*:\\s*\\S+", ln))
    if (length(m) > 0) declared = sub("^#\\s*reference\\s*:\\s*", "", m)
  }
  list(n_comment = n_comment, reference = declared)
}

# Load a gene panel TSV: `gene` column required; chrom/start/end all-or-nothing; a coordinate panel declaring another reference errors, one declaring none loads as ""
load_gene_panel = function(path, reference = NULL) {
  if (!file.exists(path)) stop("Gene panel file not found: ", path)

  hdr = .panel_header(path)
  dt = tryCatch(
    fread(path, header = TRUE, sep = "\t", fill = TRUE, skip = hdr$n_comment),
    error = function(e) fread(path, header = FALSE, sep = "\t", fill = TRUE,
                              skip = hdr$n_comment)
  )
  if (nrow(dt) == 0 && ncol(dt) == 0) stop("Gene panel file is empty: ", path)

  # copy(): setnames() rewrites the names vector in place
  orig_names = copy(names(dt))
  setnames(dt, tolower(names(dt)))
  setnames(dt, old = c("chr", "chromosome"), new = c("chrom", "chrom"), skip_absent = TRUE)

  if ("gene" %in% names(dt)) {
    genes = as.character(dt[["gene"]])
  } else if (ncol(dt) == 1) {
    # Headerless one-column list: fread consumed the first symbol as the column name, put it back
    genes = c(orig_names[1], as.character(dt[[1]]))
  } else {
    genes = as.character(dt[[1]])
  }
  genes = trimws(genes)
  keep  = nzchar(genes) & !is.na(genes) & genes != "-"

  coord_cols = c("chrom", "start", "end")
  present    = intersect(coord_cols, names(dt))
  if (length(present) > 0 && length(present) < 3) {
    stop("Gene panel ", path, " carries coordinate column(s) ",
         paste(present, collapse = ", "), " but not all of ",
         paste(coord_cols, collapse = ", "),
         ". Supply all three, or none for symbol-only matching.")
  }
  has_coords = length(present) == 3

  # A `reference` column is an alternative to the "# reference:" comment line.
  declared = hdr$reference
  if (is.na(declared) && "reference" %in% names(dt)) {
    vals = unique(trimws(as.character(dt[["reference"]])))
    vals = vals[nzchar(vals) & !is.na(vals)]
    if (length(vals) > 1)
      stop("Gene panel ", path, " declares more than one reference: ",
           paste(vals, collapse = ", "))
    if (length(vals) == 1) declared = vals
  }
  declared_norm = normalise_reference_name(declared)
  want          = normalise_reference_name(reference)

  # Only coordinate panels are reference-specific; a symbol-only panel is agnostic.
  if (has_coords && !is.na(declared_norm) && !is.na(want) && declared_norm != want) {
    stop("Gene panel ", path, " declares reference '", declared_norm,
         "' but the report is being rendered against '", want,
         "'. Panel coordinates are only valid for the reference they were built on.")
  }

  out = list(
    name       = .split_panel_filename(path)$name,
    path       = path,
    reference  = if (is.na(declared_norm)) "" else declared_norm,
    has_coords = has_coords,
    genes      = unique(genes[keep])
  )

  if (has_coords) {
    chrom = ensure_chr_prefix(trimws(as.character(dt[["chrom"]])))
    start = suppressWarnings(as.integer(dt[["start"]]))
    end   = suppressWarnings(as.integer(dt[["end"]]))
    bad   = keep & (is.na(chrom) | !nzchar(chrom) | is.na(start) | is.na(end))
    if (any(bad))
      stop("Gene panel ", path, " has missing or non-numeric coordinates for: ",
           paste(utils::head(genes[bad], 5), collapse = ", "),
           if (sum(bad) > 5) paste0(" (and ", sum(bad) - 5L, " more)") else "")
    out$chrom = chrom[keep]
    out$start = start[keep]
    out$end   = end[keep]
    # `genes` is deduplicated; the interval vectors are not (one symbol can carry several loci)
    out$interval_gene = genes[keep]
  }

  out
}

# Panel intervals as a keyed data.table, or NULL for a symbol-only panel.
panel_intervals = function(panel) {
  if (is.null(panel) || !isTRUE(panel$has_coords)) return(NULL)
  dt = data.table(gene  = as.character(panel$interval_gene),
                  chrom = as.character(panel$chrom),
                  start = as.integer(panel$start),
                  end   = as.integer(panel$end))
  setkey(dt, chrom, start, end)
  dt
}

# Is a --gene-panel argument the "no filtering" sentinel?
is_no_gene_panel = function(panel_arg) {
  is.null(panel_arg) || length(panel_arg) != 1 || is.na(panel_arg) ||
    identical(tolower(trimws(panel_arg)), "none")
}

# Path of a builtin panel, preferring the `reference`-specific variant; NULL if none
builtin_panel_path = function(assets_dir, name, reference = NULL) {
  ref = normalise_reference_name(reference)
  dir = file.path(assets_dir, "gene_lists")
  candidates = c(if (!is.na(ref)) file.path(dir, paste0(name, ".", ref, ".tsv")),
                 file.path(dir, paste0(name, ".tsv")))
  hit = candidates[file.exists(candidates)]
  if (length(hit) > 0) hit[1] else NULL
}

# Resolve a --gene-panel arg: "none", a builtin name, or a TSV path; anything else errors
resolve_gene_panel = function(panel_arg, assets_dir, reference = NULL) {
  if (is_no_gene_panel(panel_arg)) return(NULL)
  builtin = builtin_panel_path(assets_dir, panel_arg, reference)
  if (!is.null(builtin)) return(load_gene_panel(builtin, reference))
  if (file.exists(panel_arg)) return(load_gene_panel(panel_arg, reference))
  stop("Gene panel not found (tried builtin '", panel_arg, "' and as file path)")
}

# Load all builtin panels, resolving reference-specific files to one entry and skipping panels not shipped for this reference
load_all_gene_panels = function(assets_dir, reference = NULL) {
  tsv_files = Sys.glob(file.path(assets_dir, "gene_lists", "*.tsv"))
  if (length(tsv_files) == 0) return(list())

  meta = lapply(tsv_files, .split_panel_filename)
  ref  = normalise_reference_name(reference)

  panels = list()
  for (nm in unique(vapply(meta, `[[`, character(1), "name"))) {
    idx  = which(vapply(meta, `[[`, character(1), "name") == nm)
    refs = vapply(meta[idx], function(m) m$reference, character(1))
    pick = if (!is.na(ref) && any(refs == ref, na.rm = TRUE)) idx[which(refs == ref)[1]]
           else if (any(is.na(refs))) idx[which(is.na(refs))[1]]
           else NA_integer_
    if (is.na(pick)) {
      message("Gene panel '", nm, "' ships only for reference(s) ",
              paste(unique(refs), collapse = ", "), " — not offered for ",
              if (is.na(ref)) "an unknown reference" else ref)
      next
    }
    p = tryCatch(load_gene_panel(tsv_files[pick], reference),
      error = function(e) {
        message("Skipping gene panel '", nm, "': ", conditionMessage(e)); NULL })
    if (!is.null(p)) panels[[nm]] = p
  }
  panels
}

# Non-colliding key for a user panel: the first collision gets "-custom", further ones are numbered
unique_panel_name = function(nm, taken) {
  if (!(nm %in% taken)) return(nm)
  cand = paste0(nm, "-custom")
  i = 1L
  while (cand %in% taken) {
    i = i + 1L
    cand = paste0(nm, "-custom", i)
  }
  cand
}

# Resolve selected panel keys to panel objects, re-read from disk (the YAML round trip list-ifies the vectors); not tryCatch-wrapped so an unresolvable panel fails the render
resolve_selected_panels = function(keys, all_panels, assets_dir, reference = NULL) {
  keys = setdiff(as.character(unlist(keys)), "__all__")
  keys = keys[!is.na(keys) & nzchar(keys)]
  if (length(keys) == 0) return(list())
  out = list()
  for (k in unique(keys)) {
    p    = if (!is.null(all_panels)) all_panels[[k]] else NULL
    path = if (!is.null(p) && !is.null(p$path)) as.character(p$path)[1] else NA_character_
    out[[k]] = if (!is.na(path) && file.exists(path)) load_gene_panel(path, reference)
               else resolve_gene_panel(k, assets_dir, reference)
  }
  out[!vapply(out, is.null, logical(1))]
}

# Pull every occurrence of a repeatable flag out of argv (optparse has no action="append"); accepts `--flag value` and `--flag=value`
extract_repeated_option = function(args, flag) {
  args = as.character(args)
  vals = character(0)
  rest = character(0)
  i = 1L
  while (i <= length(args)) {
    a = args[i]
    if (identical(a, flag)) {
      if (i == length(args)) stop(flag, " requires a value")
      vals = c(vals, args[i + 1L])
      i = i + 2L
    } else if (startsWith(a, paste0(flag, "="))) {
      vals = c(vals, substring(a, nchar(flag) + 2L))
      i = i + 1L
    } else {
      rest = c(rest, a)
      i = i + 1L
    }
  }
  list(values = vals, rest = rest)
}

# ---- Small JS serialisation helpers: the R table and the JS indexing it are generated from the same object ----

js_quote = function(x) paste0('"', gsub('"', '\\\\"', as.character(x)), '"')

# {"col":0,"other":1} — column name to zero-based index, for a rownames=FALSE DT.
js_col_index_map = function(nms) {
  if (length(nms) == 0) return("{}")
  paste0("{", paste0(js_quote(nms), ":", seq_along(nms) - 1L, collapse = ","), "}")
}

# One JS number literal per element; element-wise so format() cannot pad or switch to scientific notation
js_num = function(x) {
  vapply(x, function(v) {
    if (is.na(v)) "null" else format(v, scientific = FALSE, trim = TRUE)
  }, character(1), USE.NAMES = FALSE)
}

.js_cell = function(v) if (is.numeric(v)) js_num(v) else ifelse(is.na(v), "null", js_quote(v))

# JS array literal from an atomic vector (numbers bare, NA as null); hand-written to avoid a jsonlite dependency
js_vec = function(x) {
  if (length(x) == 0) return("[]")
  paste0("[", paste(.js_cell(x), collapse = ","), "]")
}

# JS array of arrays, one per row of `dt`, columns in `cols` order (row-major to keep the payload small)
js_rows = function(dt, cols) {
  if (is.null(dt) || nrow(dt) == 0) return("[]")
  cells = lapply(cols, function(cl) .js_cell(dt[[cl]]))
  rows = do.call(paste, c(cells, sep = ","))
  paste0("[[", paste(rows, collapse = "],["), "]]")
}

# VEP impact severity order, as used by the styleEqual() palettes
IMPACT_LEVELS = c("HIGH", "MODERATE", "LOW", "MODIFIER")

# Distinct-value bounds for offering a dropdown; outside them a column keeps its text filter
FACET_MIN_VALUES = 2L
FACET_MAX_VALUES = 200L

# Distinct values with row counts for the tickbox column filters (assets/js/facet_filter.js): `seps` splits a column into tokens, `levels` fixes value order; returns a JS object literal {col: {sep, values: [[value, n], ...]}} with null as the "(none)" bucket. Columns absent or outside the FACET_*_VALUES bounds are omitted with a message
js_facet_defs = function(dt, cols, seps = character(0), levels = list()) {
  if (is.null(dt) || nrow(dt) == 0 || length(cols) == 0) return("{}")
  entries = character(0)

  for (nm in cols) {
    if (!nm %in% names(dt)) {
      message("Facet column '", nm, "' is not in the table - no value filter for it.")
      next
    }
    sp = if (nm %in% names(seps)) seps[[nm]] else NA_character_
    v  = as.character(dt[[nm]])

    if (!is.na(sp) && nzchar(sp)) {
      lst  = strsplit(v, sp, fixed = TRUE)
      lens = lengths(lst)
      # A cell that splits into nothing at all still contributes its row to the NA bucket.
      lst[lens == 0L] = NA_character_
      d = data.table(row = rep.int(seq_along(v), pmax(lens, 1L)),
                     tok = trimws(unlist(lst, use.names = FALSE)))
      d = unique(d, by = c("row", "tok"))   # "a,a" counts its row once
    } else {
      d = data.table(row = seq_along(v), tok = trimws(v))
    }
    d[is.na(tok) | !nzchar(tok), tok := NA_character_]

    cnt   = d[, .N, by = tok]
    n_val = nrow(cnt[!is.na(tok)])
    if (n_val < FACET_MIN_VALUES || n_val > FACET_MAX_VALUES) {
      message("Facet column '", nm, "': ", n_val,
              " distinct value(s) - keeping the plain text filter.")
      next
    }

    # Fixed levels first, then descending count, then alphabetical; the NA bucket sorts last
    lv   = intersect(as.character(levels[[nm]]), cnt$tok[!is.na(cnt$tok)])
    rank = ifelse(is.na(cnt$tok), length(lv) + 2L,
                  ifelse(cnt$tok %in% lv, match(cnt$tok, lv), length(lv) + 1L))
    cnt  = cnt[order(rank, -N, tok)]

    vals = sprintf("[%s,%s]",
                   ifelse(is.na(cnt$tok), "null", js_quote(cnt$tok)),
                   js_num(cnt$N))
    entries = c(entries,
                sprintf('%s:{"sep":%s,"values":[%s]}',
                        js_quote(nm),
                        if (is.na(sp) || !nzchar(sp)) "null" else js_quote(sp),
                        paste(vals, collapse = ",")))
  }

  paste0("{", paste(entries, collapse = ","), "}")
}

# Format a number for human-readable display
fmt_bp = function(x) {
  x = as.numeric(x)
  ifelse(abs(x) >= 1e6, paste0(round(x / 1e6, 1), " Mb"),
    ifelse(abs(x) >= 1e3, paste0(round(x / 1e3, 1), " kb"),
      paste0(x, " bp")))
}

# Embed a PNG as a framed figure (.report-figure in report.scss) with an optional caption
embed_png = function(path, max_width = "900px", caption = NULL) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  b64 = base64enc::base64encode(path)
  htmltools::tags$figure(
    class = "report-figure",
    htmltools::tags$img(
      src   = paste0("data:image/png;base64,", b64),
      style = paste0("max-width:", max_width, "; display:block; margin:auto;")
    ),
    if (!is.null(caption)) htmltools::tags$figcaption(class = "report-figure__caption", caption)
  )
}

# data: URI for a Wakhan Plotly HTML file, with a resize script injected so the plot fills the iframe
wakhan_plot_datauri = function(path) {
  html = paste(readLines(path, warn = FALSE), collapse = "\n")
  # Clear Plotly's inline fixed div size before relayout({autosize:true}) + Plots.resize()
  resize_script = "
<script>
window.addEventListener('load', function () {
  function resizeAll() {
    if (!window.Plotly) return;
    document.querySelectorAll('.plotly-graph-div').forEach(function (gd) {
      gd.style.width = '100%';
      gd.style.height = '100%';
      Plotly.relayout(gd, {autosize: true});
      Plotly.Plots.resize(gd);
    });
  }
  resizeAll();
  window.addEventListener('resize', resizeAll);
});
</script>
"
  if (grepl("</body>", html, fixed = TRUE)) {
    html = sub("</body>", paste0(resize_script, "</body>"), html, fixed = TRUE)
  } else {
    html = paste0(html, resize_script)
  }
  paste0("data:text/html;base64,", base64enc::base64encode(charToRaw(html)))
}

# Embed a self-contained HTML file (e.g. a standalone Plotly plot) as an inline iframe.
embed_html_iframe = function(path, height = "780px") {
  if (is.null(path) || !file.exists(path)) return(NULL)
  htmltools::tags$iframe(
    class = "report-figure__frame",
    src   = wakhan_plot_datauri(path),
    style = paste0("width:100%; height:", height, "; border:none;")
  )
}

# Wakhan ranked CN plots as a self-contained tab widget (Quarto's panel-tabset breaks on raw iframe HTML before a heading); panes after the first load lazily via data-src
render_wakhan_cn_tabs = function(plots) {
  if (length(plots) == 0) return(NULL)

  ids = paste0("wakhan-cn-pane-", seq_along(plots))

  buttons = lapply(seq_along(plots), function(i) {
    p = plots[[i]]
    htmltools::tags$button(
      class = if (i == 1) "wakhan-cn-tab active" else "wakhan-cn-tab",
      `data-target` = ids[i],
      paste0("Rank ", p$rank, " — purity ", p$purity, ", ploidy ", p$ploidy)
    )
  })

  panes = lapply(seq_along(plots), function(i) {
    p = plots[[i]]
    uri = wakhan_plot_datauri(p$plot)
    iframe = if (i == 1) {
      htmltools::tags$iframe(class = "report-figure__frame", src = uri,
                             style = "width:100%; height:780px; border:none;")
    } else {
      htmltools::tags$iframe(class = "report-figure__frame", `data-src` = uri,
                             style = "width:100%; height:780px; border:none;")
    }
    htmltools::tags$div(
      class = if (i == 1) "wakhan-cn-pane active" else "wakhan-cn-pane",
      id    = ids[i],
      iframe
    )
  })

  # Styled by the .wakhan-cn-tab* rules in report.scss
  htmltools::tagList(
    htmltools::tags$div(class = "wakhan-cn-tabs__nav", buttons),
    htmltools::tags$div(class = "wakhan-cn-tabs__panes", panes),
    htmltools::tags$script(htmltools::HTML("
      document.querySelectorAll('.wakhan-cn-tabs__nav').forEach(function (nav) {
        nav.querySelectorAll('.wakhan-cn-tab').forEach(function (btn) {
          btn.addEventListener('click', function () {
            const container = nav.nextElementSibling;
            nav.querySelectorAll('.wakhan-cn-tab').forEach(function (b) { b.classList.remove('active'); });
            btn.classList.add('active');
            container.querySelectorAll('.wakhan-cn-pane').forEach(function (p) { p.classList.remove('active'); });
            const pane = document.getElementById(btn.dataset.target);
            pane.classList.add('active');
            const iframe = pane.querySelector('iframe[data-src]');
            if (iframe) {
              iframe.src = iframe.dataset.src;
              iframe.removeAttribute('data-src');
            }
          });
        });
      });
    "))
  )
}

# Coding TMB from a variant_table; consequence may be comma-joined; denominator_mb defaults to 30
compute_tmb = function(variant_table, denominator_mb = 30) {
  nonsyn_terms = c(
    "missense_variant", "frameshift_variant", "stop_gained", "stop_lost",
    "start_lost", "inframe_insertion", "inframe_deletion",
    "splice_acceptor_variant", "splice_donor_variant", "protein_altering_variant"
  )
  if (is.null(variant_table) || nrow(variant_table) == 0) {
    return(list(n_nonsyn = NA_integer_, tmb = NA_real_, denominator_mb = denominator_mb))
  }
  is_nonsyn = vapply(variant_table$consequence, function(csq) {
    if (is.na(csq) || csq == "") return(FALSE)
    any(trimws(unlist(strsplit(csq, ","))) %in% nonsyn_terms)
  }, logical(1))
  n_nonsyn = sum(is_nonsyn, na.rm = TRUE)
  list(
    n_nonsyn       = n_nonsyn,
    tmb            = round(n_nonsyn / denominator_mb, 2),
    denominator_mb = denominator_mb
  )
}
