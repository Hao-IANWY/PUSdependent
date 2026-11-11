#!/usr/bin/env bash
# data.sh — collapse WT per-tissue tRNA positive sites into a representative
#           nuclear tRNA x canonical 0-76 position ratio table for the
#           plot.R heatmap.
#
# Design:
#   1. Representative tRNAs come from
#      tRNA/count/ncRNA/family_top/family_top_tRNA.tsv (the highest mean_CPM
#      member of each anticodon family, i.e. the representative by
#      expression).
#   2. Sprinzl coordinates are not converted here: the Sprinzl_pos column of
#      the pre-annotated WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv
#      under tRNA/STAT/fisher/ncRNA is read directly. Each copy uses its own
#      segment annotation, so copies whose length differs from the family
#      representative are not excluded.
#   3. Rows stay at the family level (all copies of a family fold into one
#      row); the representative only supplies the row's Representative_ID
#      label.
#
# Inputs:
#   tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv
#       chr, pos, Sprinzl_pos, strand, ratio, depth
#   tRNA/count/ncRNA/family_top/family_top_tRNA.tsv
#       feature_id, family, family_size, mean_CPM, ...
#
# Outputs:
#   representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv
#       tissue, Rep_ID, Representative_ID, chr, pos, strand, ratio, depth, canonical_pos
#   trna_landscape_summary.tsv
#       tissue, sites_in, excluded_nonNu, unmapped, sites_mapped, reps_covered, cells
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
readonly INPUT_DIR="$REPO_ROOT/tRNA/STAT/fisher/ncRNA"
readonly REPS="$REPO_ROOT/tRNA/count/ncRNA/family_top/family_top_tRNA.tsv"
readonly COLUMNS=$'tissue\tRep_ID\tRepresentative_ID\tchr\tpos\tstrand\tratio\tdepth\tcanonical_pos'

log() { printf '[data_prep] %s\n' "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

[[ -s "$REPS" ]] || error "Missing family_top table: $REPS"

readonly SUMMARY="$SCRIPT_DIR/trna_landscape_summary.tsv"
readonly COUNTS="$SCRIPT_DIR/.counts.tmp"
trap 'rm -f "$COUNTS"' EXIT

printf 'tissue\tsites_in\texcluded_nonNu\tunmapped\tsites_mapped\treps_covered\tcells\n' > "$SUMMARY"

for tissue in Root Leaf Flower; do
    input="$INPUT_DIR/WT_${tissue}_sites.positive.Sprinzl.tsv"
    [[ -s "$input" ]] || error "Missing annotated input: $input (run add_sprinzl_pos.py first)"
    prefix="$(printf '%s' "$tissue" | tr '[:upper:]' '[:lower:]')"
    out="$SCRIPT_DIR/representative_nuclear_trna_mapped_sites_${prefix}.tsv"

    printf '%s\n' "$COLUMNS" > "$out"
    awk -F'\t' -v OFS='\t' -v tissue="$tissue" -v reps="$REPS" -v counts="$COUNTS" '
    # Family representative table: nuclear tRNAs only; the family name is the
    # Rep_ID (e.g. Nu-Arg-TCG)
    FILENAME == reps {
        if (FNR > 1 && $1 ~ /^tRNA-Nu-/) REP[$2] = $1
        next
    }
    FNR == 1 { next }
    {
        sites_in++
        chr = $1
        if (chr !~ /^tRNA-Nu-/) { excluded++; next }        # nuclear tRNAs only
        fam = chr
        sub(/^tRNA-Nu-/, "Nu-", fam)
        sub(/-[0-9]+$/, "", fam)
        if (!(fam in REP)) { excluded++; next }             # families absent from family_top

        canon = $3
        if (canon == "" || canon !~ /^[0-9]+$/ || canon + 0 < 0 || canon + 0 > 76) {
            unmapped++; next                                # Sprinzl_pos missing or out of range
        }

        print tissue, fam, REP[fam], chr, $2, $4, sprintf("%.10g", $5), sprintf("%.10g", $6), canon + 0
        mapped++
        seen_rep[fam] = 1
        seen_cell[fam SUBSEP canon] = 1
    }
    END {
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", tissue, sites_in + 0, excluded + 0,
               unmapped + 0, mapped + 0, length(seen_rep), length(seen_cell) > counts
    }
    ' "$REPS" "$input" | sort -t$'\t' -k2,2 -k9,9n -k5,5n >> "$out"

    counts_line="$(cat "$COUNTS")"
    printf '%s\n' "$counts_line" >> "$SUMMARY"
    log "$(printf '%s' "$counts_line" | tr '\t' ' ') -> $(basename "$out")"
done

log "Wrote: $(basename "$SUMMARY")"
log "Now drawing figures"
Rscript "$SCRIPT_DIR/plot.R"
