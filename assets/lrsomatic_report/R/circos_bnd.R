suppressPackageStartupMessages({
  library(data.table)
})

# Breakend-only circos cross-linked to the SV table; R selects the data, assets/js/bnd_circos.js draws it so sectors re-lay-out per filter

# Payload cap: real samples carry 50-115 arcs, an unfiltered run can carry ~13k
BND_CIRCOS_MAX_LINKS = 2000L

# BND classes with two loci to draw between (a single breakend has no partner)
BND_CIRCOS_CLASSES = c("translocation", "intra-chr breakend")

# ---- Data selection ------------------------------------------------------

# Drawable rearrangements: BNDs with both loci on plotted chromosomes; feeds both bnd_panel_genes() and bnd_circos_data()
bnd_links = function(sv_table, chromosomes) {
  cols = c("id", "svclass", "chrom_a", "pos_a", "chrom_b", "pos_b")
  empty = data.table(id = character(), svclass = character(),
                     chrom_a = character(), pos_a = integer(),
                     chrom_b = character(), pos_b = integer())
  if (is.null(sv_table) || nrow(sv_table) == 0 || !all(cols %in% names(sv_table)))
    return(empty)

  b = as.data.table(sv_table)[, ..cols]
  b = b[svclass %in% BND_CIRCOS_CLASSES &
        !is.na(chrom_a) & !is.na(pos_a) & !is.na(chrom_b) & !is.na(pos_b) &
        chrom_a %in% chromosomes & chrom_b %in% chromosomes]
  if (nrow(b) == 0) return(empty)

  # Deterministic order, so the payload is byte-identical across renders of the same data.
  setorder(b, chrom_a, pos_a, chrom_b, pos_b, id)
  b[]
}

# Panel genes within `window` of a drawn breakend, using the same test as sv_panel_hits() so the two cannot drift
bnd_panel_genes = function(links, all_panels, window = SV_PANEL_WINDOW_BND) {
  empty = data.table(chrom = character(), start = integer(), end = integer(),
                     gene = character(), panels = character())
  if (is.null(links) || nrow(links) == 0) return(empty)
  if (is.null(all_panels) || length(all_panels) == 0) return(empty)

  q = rbindlist(list(
    data.table(chrom = links$chrom_a, pos = as.numeric(links$pos_a)),
    data.table(chrom = links$chrom_b, pos = as.numeric(links$pos_b))
  ))
  q = unique(q[!is.na(chrom) & !is.na(pos)])
  if (nrow(q) == 0) return(empty)
  q[, `:=`(start = pmax(pos - window, 0), end = pos + window)]

  out = list()
  for (nm in names(all_panels)) {
    iv = panel_intervals(all_panels[[nm]])
    if (is.null(iv) || nrow(iv) == 0) next
    iv = copy(iv)[, `:=`(start = as.numeric(start), end = as.numeric(end))]
    setkey(iv, chrom, start, end)
    ov = data.table::foverlaps(q[, .(chrom, start, end)], iv,
                               by.x = c("chrom", "start", "end"),
                               type = "any", nomatch = NULL)
    if (nrow(ov) == 0) next
    out[[nm]] = unique(ov[, .(chrom, start, end, gene, panel = nm)])
  }
  if (length(out) == 0) return(empty)

  g = rbindlist(out)
  g = g[, .(panels = paste(sort(unique(panel)), collapse = ",")),
        by = .(chrom, start, end, gene)]
  setorder(g, chrom, start, gene)
  g[, `:=`(start = as.integer(start), end = as.integer(end))]
  g[, .(chrom, start, end, gene, panels)]
}


# ---- Payload for the client-side plot -------------------------------------

# Payload for the client-side breakend circos: list(data, n_links, n_genes, chroms, reason); data is NULL when there is nothing, or too much, to draw
bnd_circos_data = function(links, genes = NULL, cytobands, chrom_lengths, chromosomes) {

  fail = function(reason) list(data = NULL, n_links = 0L, n_genes = 0L,
                               chroms = character(0), reason = reason)

  if (is.null(links) || nrow(links) == 0)
    return(fail("No breakends with two mapped loci — nothing to draw arcs between."))
  if (nrow(links) > BND_CIRCOS_MAX_LINKS)
    return(fail(sprintf(paste("%s breakend arcs is more than this plot can show",
                              "(limit %s). Select a gene panel above for a filtered",
                              "report, or read the table below."),
                        nrow(links), BND_CIRCOS_MAX_LINKS)))

  b = copy(as.data.table(links))
  chroms = chromosomes[chromosomes %in% unique(c(b$chrom_a, b$chrom_b))]
  chroms = chroms[chroms %in% unique(cytobands$chrom)]
  b = b[chrom_a %in% chroms & chrom_b %in% chroms]
  if (nrow(b) == 0 || length(chroms) == 0)
    return(fail("No breakends on the chromosomes this report plots."))
  chroms = chromosomes[chromosomes %in% unique(c(b$chrom_a, b$chrom_b))]

  lens = chrom_lengths[chroms]
  # Clamp a locus past its sector end rather than dropping the row
  b[, `:=`(pos_a = pmin(pmax(as.numeric(pos_a), 1), lens[chrom_a]),
           pos_b = pmin(pmax(as.numeric(pos_b), 1), lens[chrom_b]))]

  cy = as.data.table(cytobands)[chrom %in% chroms]
  setorder(cy, chrom, start)

  g = if (is.null(genes)) data.table() else as.data.table(genes)[chrom %in% chroms]
  if (nrow(g) > 0) setorder(g, chrom, start, gene)

  list(
    data = list(
      chromosomes = js_vec(chroms),
      lengths     = js_vec(as.numeric(lens)),
      cytobands   = js_rows(cy, c("chrom", "start", "end", "stain")),
      links       = js_rows(b,  c("id", "svclass", "chrom_a", "pos_a", "chrom_b", "pos_b")),
      genes       = if (nrow(g) > 0) js_rows(g, c("chrom", "start", "end", "gene", "panels"))
                    else "[]"
    ),
    n_links = nrow(b),
    n_genes = nrow(g),
    chroms  = chroms,
    reason  = NULL
  )
}

# The <script> that hands the payload to the client, as one window.BND_DATA literal.
bnd_circos_script = function(res) {
  if (is.null(res$data)) return("")
  d = res$data
  paste0(
    "<script>window.BND_DATA={",
    "chromosomes:", d$chromosomes, ",",
    "lengths:",     d$lengths,     ",",
    "cytobands:",   d$cytobands,   ",",
    "links:",       d$links,       ",",
    "genes:",       d$genes,
    "};</script>")
}
