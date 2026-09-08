suppressPackageStartupMessages({
  library(data.table)
})

# Load cytobands for a given reference; returns data.frame suitable for circlize
load_cytobands = function(reference, assets_dir) {
  ref = tolower(reference)
  path = file.path(assets_dir, "references", ref, "cytobands.tsv")
  if (!file.exists(path)) stop("No cytobands for reference '", ref, "': ", path)
  dt = fread(path, header = FALSE, sep = "\t",
             col.names = c("chrom", "start", "end", "name", "stain"))
  as.data.frame(dt)
}

# Load chromosome lengths; returns named integer vector (name = chrom, value = length)
load_chrom_lengths = function(reference, assets_dir) {
  ref = tolower(reference)
  path = file.path(assets_dir, "references", ref, "chrom_lengths.tsv")
  if (!file.exists(path)) stop("No chrom_lengths for reference '", ref, "': ", path)
  dt = fread(path, header = FALSE, sep = "\t", col.names = c("chrom", "length"))
  setNames(as.integer(dt$length), dt$chrom)
}

# Auto-detect reference from VCF/VEP headers (chr1 length: CHM13v2 248387328, GRCh38 248956422; assembly/genome_build lines)
detect_reference = function(vcf_file) {
  if (!file.exists(vcf_file)) {
    message("Cannot auto-detect reference: file not found, defaulting to t2t")
    return("t2t")
  }
  con = gzfile(vcf_file, "rb")
  on.exit(close(con))
  header_lines = character(0)
  for (i in seq_len(2000)) {
    line = tryCatch(readLines(con, n = 1, warn = FALSE), error = function(e) character(0))
    if (length(line) == 0 || !startsWith(line, "##")) break
    header_lines = c(header_lines, line)
  }

  # 1. Check VEP "## assembly version" line
  asm_line = grep("assembly version|genome_build|assembly=", header_lines,
                  value = TRUE, ignore.case = TRUE)
  if (length(asm_line) > 0) {
    asm = tolower(paste(asm_line, collapse = " "))
    if (grepl("t2t|chm13", asm))   return("t2t")
    if (grepl("grch38|hg38|38", asm)) return("hg38")
  }

  # 2. Check ##contig chr1 length (standard VCF)
  contig_chr1 = grep("ID=chr1[^0-9].*length=|ID=1[^0-9].*length=",
                     header_lines, value = TRUE, perl = TRUE)
  if (length(contig_chr1) > 0) {
    len = as.integer(sub(".*length=([0-9]+).*", "\\1", contig_chr1[1]))
    if (!is.na(len)) {
      if (abs(len - 248387328L) < 1000L) return("t2t")
      if (abs(len - 248956422L) < 1000L) return("hg38")
    }
  }

  message("Could not determine reference from file headers, defaulting to t2t")
  "t2t"
}

# Build the chromosome list for plotting based on sex
chromosomes_for_sex = function(sex) {
  sex = tolower(trimws(sex))
  autosomes = paste0("chr", 1:22)
  if (sex %in% c("male", "xy")) c(autosomes, "chrX", "chrY")
  else c(autosomes, "chrX")
}
