#!/usr/bin/env bash
#
# Reshape the REVEL and EVE releases into the form VEP expects. The pipeline runs
# this itself; run it by hand to pre-build files for --vep_revel / --vep_eve.
#
# Usage:
#   bin/prepare_vep_plugin_data.sh revel <zip-or-unpacked-dir> <outdir>
#   bin/prepare_vep_plugin_data.sh eve <eve-vcf-dir> <outdir>
#
# Requires: bgzip and tabix (htslib), awk, sort, and unzip for the REVEL zip.

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

need() {
    for tool in "$@"; do
        command -v "$tool" >/dev/null || die "$tool not found on PATH"
    done
}

usage() {
    sed -n '2,10p' "$0" | sed 's/^#\{1,2\} \{0,1\}//'
    exit "${1:-1}"
}

# REVEL ships comma-separated, sorted on its GRCh37 column (2); GRCh38 (column 3) needs a re-sort and its own index.
prep_revel() {
    local src=${1:-}
    local outdir=${2:-.}
    [[ -n $src ]] || die "usage: $0 revel <revel-v1.3_all_chromosomes.zip|unpacked-dir> [outdir]"
    [[ -r $src ]] || die "cannot read $src"
    need bgzip tabix awk sort

    mkdir -p "$outdir"
    local work
    work=$(mktemp -d "${TMPDIR:-/tmp}/revel.XXXXXX")
    trap 'rm -rf "$work"' RETURN

    # Either the release zip or a directory it has already been unpacked into
    local searchdir=$src
    if [[ -f $src ]]; then
        need unzip
        echo ">> unpacking $src" >&2
        unzip -o -q -d "$work" "$src"
        searchdir=$work
    elif [[ ! -d $src ]]; then
        die "$src is neither a zip file nor a directory"
    fi

    # -L because the pipeline passes a Nextflow-staged directory, which is a symlink;
    # plain find reports it and never descends.
    local raw
    raw=$(find -L "$searchdir" -name 'revel_with_transcript_ids' -o -name 'revel_all_chromosomes.csv' | head -n1)
    [[ -n $raw ]] || die "could not find the REVEL table in $src"
    echo "   using $(basename "$raw")" >&2

    local out="$outdir/revel_grch38.tsv.gz"
    echo ">> converting to tab-separated and re-sorting on the GRCh38 column" >&2
    {
        # header row, tab-separated and '#'-prefixed so tabix treats it as a comment
        head -n1 "$raw" | tr ',' '\t' | sed '1s/^/#/'
        # data rows: drop those with no GRCh38 position, then sort on it
        tail -n +2 "$raw" | tr ',' '\t' | awk -F'\t' '$3 != "." && $3 != ""' \
          | sort -T "$work" -k1,1 -k3,3n
    } | bgzip -c > "$out"

    tabix -f -s 1 -b 3 -e 3 -c '#' "$out"
    echo ">> wrote $out and $out.tbi" >&2
    echo >&2
    echo "Pass these to the pipeline with:" >&2
    echo "  --vep_revel $out --vep_revel_tbi $out.tbi" >&2
}

# EVE ships one VCF per protein; the plugin expects a single sorted, indexed VCF.
prep_eve() {
    local vcfdir=${1:-}
    local outdir=${2:-.}
    [[ -n $vcfdir ]] || die "usage: $0 eve <dir-of-per-protein-vcfs> [outdir]"
    [[ -d $vcfdir ]] || die "$vcfdir is not a directory"
    need bgzip tabix awk sort

    mkdir -p "$outdir"
    local work
    work=$(mktemp -d "${TMPDIR:-/tmp}/eve.XXXXXX")
    trap 'rm -rf "$work"' RETURN

    # The bulk zip nests the per-protein files under vcf_files_missense_mutations.
    local src=$vcfdir
    if [[ -d "$vcfdir/vcf_files_missense_mutations" ]]; then
        src="$vcfdir/vcf_files_missense_mutations"
    fi

    # -L throughout because the pipeline passes a Nextflow-staged directory, which is a
    # symlink; plain find reports it and never descends.
    local count
    count=$(find -L "$src" -name '*.vcf' | wc -l)
    [[ $count -gt 0 ]] || die "no .vcf files found under $src"
    echo ">> merging $count per-protein VCFs from $src" >&2

    local first
    first=$(find -L "$src" -name '*.vcf' | sort | head -n1)

    local out="$outdir/eve_merged.vcf.gz"
    {
        grep '^#' "$first"
        find -L "$src" -name '*.vcf' -exec grep -hv '^#' {} + \
          | sort -T "$work" -k1,1V -k2,2n
    } | bgzip -c > "$out"

    tabix -f -p vcf "$out"
    echo ">> wrote $out and $out.tbi" >&2
    echo >&2
    echo "Pass these to the pipeline with:" >&2
    echo "  --vep_eve $out --vep_eve_tbi $out.tbi" >&2
}

case "${1:-}" in
    revel)         shift; prep_revel "$@" ;;
    eve)           shift; prep_eve "$@" ;;
    -h|--help|help) usage 0 ;;
    "")            die "no resource given. Try: $0 --help" ;;
    *)             die "unknown resource '$1'. Expected revel or eve." ;;
esac
