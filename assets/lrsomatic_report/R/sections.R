# Section-module contract: id, title, locate(sample_dir, sample_id), parse(inputs, section_data); see CLAUDE.md

SECTIONS = list()

register_section = function(descriptor) {
  SECTIONS[[descriptor$id]] <<- descriptor
}

# Quiet "nothing to show" notice; `warn` marks a genuine failure
section_notice = function(msg, warn = FALSE) {
  tags$div(class = if (warn) "section-notice section-notice--warn" else "section-notice", msg)
}

# Collapsed provenance/caveats under a table; `summary` names the contents
table_details = function(..., summary = "Details") {
  tags$details(
    class = "table-details",
    tags$summary(summary),
    tags$div(class = "table-footnote", ...)
  )
}
