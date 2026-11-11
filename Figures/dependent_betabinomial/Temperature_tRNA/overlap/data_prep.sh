#!/usr/bin/env bash
# data_prep.sh — Overlap table of "27 °C up-regulated sites" and "27 °C trua1b-dependent
#               sites" (tRNA, one table)
#
# Purpose: TSV output only, no plotting. The mRNA / rRNA version lives in the sibling
# directory Figures/dependent_betabinomial/Temperature_mRNA/overlap/
#
# Sources of the two sets:
#   27 °C up        : Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.response.tsv
#                     (direction == up_in_WT_27)
#   27 °C dependent : Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/27/trua1b/trua1b-27.dependent.tsv
#
# Site key = chr + pos + strand; the union of the two sets is flagged per site with a
# combination label (up_only / dependent_only / both) so plot.R can count the three
# classes and draw a 2-set Venn.
#
# tRNA note: the Temperature_tRNA reference is tRNA-unique_sn_snoRNA.fa, which mixes in
# sn/snoRNAs; both sets must be filtered with ^tRNA- (same criteria as
# Figures/dependent_betabinomial/tRNA).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly RESPONSE_TRNA="$ROOT_DIR/Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA"
readonly BB_TRNA="$ROOT_DIR/Temperature_tRNA/STAT/dependent_betabinomial"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

write_overlap_table() {
  local up_tsv="$1"
  local dep_tsv="$2"
  local output_tsv="$3"
  local rna_label="$4"

  awk -F'\t' -v OFS='\t' -v rna_label="$rna_label" '
    function keep() { return $1 ~ /^tRNA-/ }
    function register(k, c, p, s, ctx) {
      if (!(k in seen)) { seen[k] = 1; order[++n] = k }
      chr[k] = c; pos[k] = p; strand[k] = s
      if (ctx != "" && !(k in context)) { context[k] = ctx }
    }
    BEGIN { print "chr", "pos", "strand", "context", "site_id", "up_in_WT_27",
                  "trua1b_dependent_27", "combination" }
    ARGIND == 1 {
      if (FNR == 1 || !keep() || $9 != "up_in_WT_27") { next }
      gsub(/\r/, "")
      k = $1 SUBSEP $2 SUBSEP $3
      register(k, $1, $2, $3, $4)
      up[k] = 1
      next
    }
    ARGIND == 2 {
      if (FNR == 1 || !keep()) { next }
      gsub(/\r/, "")
      k = $1 SUBSEP $2 SUBSEP $3
      register(k, $1, $2, $3, $4)
      dep[k] = 1
      next
    }
    END {
      for (i = 1; i <= n; i++) {
        k = order[i]
        u = (k in up) ? 1 : 0
        d = (k in dep) ? 1 : 0
        combo = (u && d) ? "both" : ((u) ? "up_only" : "dependent_only")
        count[combo]++
        print chr[k], pos[k], strand[k], context[k], (chr[k] "|" pos[k] "|" strand[k]),
              u, d, combo
      }
      printf "[data_prep] %s: up_only=%d, dependent_only=%d, both=%d\n",
             rna_label, count["up_only"] + 0, count["dependent_only"] + 0,
             count["both"] + 0 > "/dev/stderr"
    }
  ' "$up_tsv" "$dep_tsv" > "$output_tsv"
}

write_pair() {
  local up_tsv="${RESPONSE_TRNA}/WT_27_vs_WT_21.response.tsv"
  local dep_tsv="${BB_TRNA}/ncRNA/27/trua1b/trua1b-27.dependent.tsv"

  require_file "$up_tsv"
  require_file "$dep_tsv"

  write_overlap_table "$up_tsv" "$dep_tsv" \
    "${SCRIPT_DIR}/tRNA_27up_vs_trua1b_dependent_overlap.tsv" "tRNA"
}

main() {
  write_pair
}

main "$@"
