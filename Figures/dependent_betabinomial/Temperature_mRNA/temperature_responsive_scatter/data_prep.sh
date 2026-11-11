#!/usr/bin/env bash
# data_prep.sh — Assemble the data tables for the WT temperature-responsive site
#                scatter plots at 21 / 27 °C (mRNA / rRNA, one each)
#
# Purpose: TSV output only, no plotting. The tRNA version lives in the sibling
# directory Figures/dependent_betabinomial/Temperature_tRNA/temperature_responsive_scatter/
#
# Data source:
#   mRNA / rRNA : Temperature/STAT/Temperature_response_betabinomial/<type>/WT_27_vs_WT_21.all-sites.tsv
#
# all-sites.tsv is the full candidate-site table (the union of Fisher-positive sites at
# both temperatures, with modification fractions and mean Mock depth recomputed per
# site); it is semantically equivalent to the on-the-fly
# temperature_responsive_scatter.tsv of the paper's temperature_responsive_sites
# workflow, so the stored product is reused here with one round of column-name
# normalization instead of recomputing from the raw fisher/adjustGap tables.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly TEMPERATURE_ALL_SITES="$ROOT_DIR/Temperature/STAT/Temperature_response_betabinomial"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# Turn all-sites.tsv into the plot table.
# Columns kept: chr pos strand context WT_21_Ratio WT_27_Ratio Delta_Ratio_27_minus_21
#               Mock_Mean_Depth_WT_21 Mock_Mean_Depth_WT_27 p_value q_value regulation
# regulation is rewritten from the upstream direction column: empty -> non-responsive,
# up_in_WT_27 -> up-regulated (higher modification at 27 °C), down_in_WT_27 ->
# down-regulated.
# Calling criteria (see each STAT/Temperature_response_betabinomial/run.sh):
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
  local rna_type
  for rna_type in mRNA rRNA; do
    local input_tsv="${TEMPERATURE_ALL_SITES}/${rna_type}/WT_27_vs_WT_21.all-sites.tsv"
    local output_tsv="${SCRIPT_DIR}/${rna_type}_temperature_responsive_scatter.tsv"
    require_file "$input_tsv"
    write_plot_table "$input_tsv" "$output_tsv"
  done
}

main "$@"
