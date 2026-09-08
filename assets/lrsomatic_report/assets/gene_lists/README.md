# Gene Panel Lists

Each file is a TSV with a required `gene` column (HGNC symbol). Coordinate columns
`chrom` (or `chr`), `start` and `end` are optional but **all-or-nothing** — a file with
some but not all three is rejected rather than quietly falling back to symbol matching.

| Panel columns | Small-variant filter | SV filter |
|---|---|---|
| `gene` only | symbol match | direct-hit symbol match on either breakend's VEP gene — no windows |
| `gene, chrom, start, end` | symbol match | coordinate match: within 1 Mb of a breakend (BND) or 100 kb of the SV span (other types) |

Coordinate matching is what makes breakend filtering reliable: whether a BND carries a
VEP gene symbol at all depends on the sample's VEP invocation (1.6%–90% of breakends
across the samples measured), so a symbol-only panel can hide the very translocations it
exists to find. Matching on coordinates needs no annotation on the row.

Optional metadata columns (`panel`, `notes`) are ignored by the loader and kept for the
reader.

## Reference declaration

Panel coordinates are only valid for the reference they were built on — matching an hg38
panel against a T2T sample produces wrong hits with no error anywhere. A
coordinate-carrying panel must therefore declare its reference, either as a leading
comment line:

```
# reference: hg38
gene	chrom	start	end
MYC	chr8	127735434	127742951
```

or as a `reference` column. A panel whose declared reference differs from the one the
report is rendered against is a **hard error**. A panel that declares none loads, but the
SV section footnote says "reference unverified". Symbol-only panels are
reference-agnostic and need no declaration.

Builtin panels ship one file per reference (`lymphoid.hg38.tsv`, `lymphoid.t2t.tsv`) and
are presented as a single selectable `lymphoid` entry, resolved against the detected
reference.

## Supplying a custom panel

```bash
Rscript bin/render_report.R \
  --sample-dir /path/to/sample \
  --sample-id MySample \
  --gene-panel /path/to/my_genes.tsv
```

A one-column file of symbols (with or without a `gene` header) is accepted, and gives
symbol-only matching. The report's "Custom…" textarea takes bare symbols, so it is
symbol-only too.

`--gene-panel` is repeatable, so several panels can be applied at once — a builtin and your
own list together, say. A variant or SV is kept if it hits any of them:

```bash
  --gene-panel lymphoid --gene-panel /path/to/my_genes.tsv
```

Each panel is registered under its filename stem; two files sharing a basename both stay
selectable, the second as `<name>-custom`. Every registered panel is a checkbox in the
report, so the reader can retick them without re-rendering.

## Bundled panels

| File | Contents |
|---|---|
| `lymphoid.hg38.tsv` | 72 recurrently mutated genes in B-cell lymphomas (DLBCL, FL, MCL, CLL, BL, MALT), GENCODE v46 gene spans |
| `lymphoid.t2t.tsv` | the same 72 genes, spans from the CHM13v2.0 RefSeq Liftoff v5.1 annotation |
| `sarcoma.hg38.tsv` | 140 soft-tissue and bone sarcoma genes — tumour suppressors, amplification targets and recurrent fusion partners — GENCODE v46 gene spans |
| `sarcoma.t2t.tsv` | the same 140 genes, spans from the CHM13v2.0 RefSeq Liftoff v5.1 annotation |

Regenerating them is mechanical — gene spans keyed on `gene_name`, taken from
`gene` features (GENCODE) or the min/max of `transcript` features (Liftoff, which has no
`gene` feature), restricted to `chr1`–`chr22`, `chrX`, `chrY`:

- hg38: `references/GRCh38.alt-masked-V2/annotation/gencode.v46.basic.annotation.gtf.gz`
- t2t: `references/chm13_v2.0_maskedY.rCRS/annotation/chm13v2.0_RefSeq_Liftoff_v5.1.gtf`

Neither annotation is keyed on current HGNC symbols throughout, so aliases in the source
gene lists were mapped by hand. For `lymphoid`: `CD20`→`MS4A1` (already present, so the
rows merged), `GEF1`→`ARHGEF1`, `HIST1H1E`→`H1-4`. For `sarcoma`, from an input list of
151 lines: `VEGFR2`→`KDR`, `VEGFR3`→`FLT4`, `MKL2`→`MRTFB`, `MGEA5`→`OGA`, plus
`HER2`→`ERBB2`, `SYT`→`SS18`, `H3F3A`→`H3-3A` and `H3F3B`→`H3-3B`, whose targets were
already listed — those merged, as did seven verbatim duplicates, leaving 140 genes.

Two `sarcoma` symbols need the T2T annotation handled specially, and both are recorded in
that file's comment header:

- `POU2AF3` — the Liftoff annotation predates the rename and carries it as `COLCA2`.
- `DUX4L10` — the Liftoff annotation has no D4Z4 paralogs at all (only `DUX4` itself), so
  this span comes from `GCF_009914755.1_T2T-CHM13v2.0_genomic.gtf.gz`, whose chromosomes
  are NCBI accessions (`NC_060925.1` = `chr1` … `NC_060948.1` = `chrY`). Its 10q26
  position agrees with the hg38 locus, so this is a genuine match rather than a guess.

A missing coordinate is a **hard error** in `load_gene_panel()`, not a per-row downgrade to
symbol matching, so a gene that resolves in one annotation and not the other has to be
either mapped or dropped from that reference's file — it cannot be left blank.
