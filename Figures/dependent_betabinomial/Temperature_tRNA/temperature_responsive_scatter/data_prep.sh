#!/usr/bin/env bash
# data_prep.sh — Assemble the data table for the WT temperature-responsive site scatter
#                plot at 21 / 27 °C (tRNA, one table)
#
# Purpose: TSV output only, no plotting. The mRNA / rRNA version lives in the sibling
# directory Figures/dependent_betabinomial/Temperature_mRNA/temperature_responsive_scatter/
#
# Data source:
#   tRNA : Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.all-sites.tsv
#
# all-sites.tsv is the full candidate-site table (the union of Fisher-positive sites at
# both temperatures, with modification fractions and mean Mock depth recomputed per site
# from adjustGap); it is semantically equivalent to the on-the-fly
# temperature_responsive_scatter.tsv of the paper's temperature_responsive_sites
# workflow, so the stored product is reused here with one round of column-name
# normalization instead of recomputing from the raw fisher/adjustGap tables.
#
# tRNA note: the Temperature_tRNA reference is tRNA-unique_sn_snoRNA.fa, which mixes in
# sn/snoRNAs; filter with ^tRNA- to keep only the tRNA fraction (same criteria as
# Figures/dependent_betabinomial/tRNA).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly TEMPERATURE_TRNA_ALL_SITES="$ROOT_DIR/Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.all-sites.tsv"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# Turn all-sites.tsv into the plot table.
# Columns kept: chr pos strand context WT_21_Ratio WT_27_Ratio Delta_Ratio_27_minus_21
#               Mock_Mean_Depth_WT_21 Mock_Mean_Depth_WT_27 p_value q_value regulation
# regulation is rewritten from the upstream direction column: empty -> non-responsive,
# up_in_WT_27 -> up-regulated (higher modification at 27 °C), down_in_WT_27 ->
# down-regulated.
# Calling criteria (see Temperature_tRNA/STAT/Temperature_response_betabinomial/run.sh):
# beta-binomial GLMM with a global theta + per-site LRT (--site-test lrt), BH q < 0.05
# (two-sided FDR), |Delta| >= 0.05, and depth >= 10 in every replicate of both
# temperature groups.
write_plot_table() {
  local input_tsv="$1"
  local output_tsv="$2"

  awk -F'\t' -v OFS='\t' '
    BEGIN { print "chr", "pos", "strand", "context", "WT_21_Ratio", "WT_27_Ratio",
                  "Delta_Ratio_27_minus_21", "Mock_Mean_Depth_WT_21", "Mock_Mean_Depth_WT_27",
                  "p_value", "q_value", "regulation" }
    NR == 1 { next }
    {
      if ($1 !~ /^tRNA-/) { next }
      regulation = "non-responsive"
      if ($13 == "up_in_WT_27") { regulation = "up-regulated" }
      else if ($13 == "down_in_WT_27") { regulation = "down-regulated" }
      print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, regulation
      n++
      if (regulation != "non-responsive") { n_responsive++ }
    }
    END { printf "[data_prep] %s: %d sites, %d responsive\n", FILENAME, n, n_responsive + 0 > "/dev/stderr" }
  ' "$input_tsv" > "$output_tsv"
}

main() {
  require_file "$TEMPERATURE_TRNA_ALL_SITES"
  write_plot_table "$TEMPERATURE_TRNA_ALL_SITES" \
    "${SCRIPT_DIR}/tRNA_temperature_responsive_scatter.tsv"
}

main "$@"
