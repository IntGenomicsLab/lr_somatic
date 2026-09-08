# Vendored: lrsomatic_report

This directory is a **vendored copy** of the standalone report tool, not a git submodule.
Do not edit it here — fix upstream, tag a release, and re-sync.

| | |
|---|---|
| Upstream | <https://github.com/ljwharbers/lrsomatic_report> |
| Release | `v1.3.2` (`75c65b2b68968f3d742dbfd1ca7ef780dbdaf770`) |
| Vendored commit | `75c65b2b68968f3d742dbfd1ca7ef780dbdaf770` (the tag itself) |
| License | MIT (see `LICENSE`) |

## Why vendored rather than a submodule

`nextflow run IntGenomicsLab/lrsomatic` clones the pipeline repository but does **not**
fetch git submodules, so a gitlink here would be an empty directory for every user who
did not hand-clone with `--recurse-submodules` — and for CI, whose checkout steps do not
pass `submodules: recursive`. Real tracked files work for both.

The tool's *dependencies* (R, Quarto and its R packages) are handled separately, by the
Wave multi-package container declared in `modules/local/lrsomaticreport/main.nf` and
built from that module's `environment.yml`.

## What is included

Only what `bin/render_report.R` needs at run time:

```
bin/  R/  templates/  assets/  LICENSE  README.md
```

Upstream `docs/`, `tests/` (both `tests/testthat/` and `tests/js/`), `recipe/` and
`CLAUDE.md` are deliberately excluded. (Upstream marks `docs/` and `tests/`
`export-ignore` in its `.gitattributes`; `recipe/` and `CLAUDE.md` are not marked, so
they must be left behind by hand when copying.)

## Re-syncing on the next upstream release

```bash
TAG=v1.4.0
git clone --depth 1 --branch "$TAG" https://github.com/ljwharbers/lrsomatic_report.git "$TMPDIR/lrr"
rm -rf assets/lrsomatic_report/{bin,R,templates,assets,LICENSE,README.md}
cp -a "$TMPDIR"/lrr/{bin,R,templates,assets,LICENSE,README.md} assets/lrsomatic_report/
```

The `rm -rf` before the copy is not optional: it is what removes files *deleted* upstream.
Copying over the top would have left the pre-v1.2.1 `assets/gene_lists/lymphoid.tsv` behind,
where `load_all_gene_panels()` would glob it as a second, reference-less `lymphoid` panel
alongside the `lymphoid.hg38.tsv` / `lymphoid.t2t.tsv` pair that replaced it.

Then update the table above, and:

- `modules/local/lrsomaticreport/environment.yml` — if upstream `recipe/meta.yaml` gained a
  dependency. The two files are kept in exact agreement; check with a `library()`/`require()`
  grep over the upstream `R/`, `bin/` and `templates/` rather than trusting the README.
- `modules/local/lrsomaticreport/main.nf` — the hard-coded version topic (the tool has no
  `--version` flag), and the container digests **only if `environment.yml` changed**.
- `modules/local/lrsomaticreport/meta.yml` — the same version string, in two places.
- `modules/local/lrsomaticreport/tests/main.nf.test.snap` — regenerate.
- `tests/{default,clair_only,deep_only,consensus,union}.nf.test.snap` — the
  `"lrsomatic_report": "<version>"` line in each.
- `CHANGELOG.md` — the vendored version named in the `#176` entry.

Nothing checks these for consistency; miss one and the snapshots go stale silently.

Rebuild the containers after any `environment.yml` change so the images and the file agree.
**Two** builds are needed and both must be updated in `main.nf`: `--singularity` produces a
Singularity-native SIF (the `oras://` reference), the default build a genuine OCI image.

```bash
wave --conda-file modules/local/lrsomaticreport/environment.yml --freeze --await
wave --conda-file modules/local/lrsomaticreport/environment.yml --freeze --await --singularity
```
