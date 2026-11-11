#!/usr/bin/env bash
# data_prep.sh — Assemble the data tables for the loss-of-modification dot plots of
#                 trua1b-dependent Ψ sites (mRNA / rRNA × 21 °C / 27 °C, 4 tables)
#
# Purpose: TSV output only, no plotting. The tRNA version lives in the sibling
# directory Figures/dependent_betabinomial/Temperature_tRNA/dependent_site_loss_effect/
#
# Data source:
#   mRNA / rRNA : Temperature/STAT/dependent_betabinomial/<type>/<temp>/trua1b/
#
# Three tables per combination:
#   trua1b_1_<temp>.dependent-sites.tsv   per-site table of allele trua1b_1 (all WT candidate sites)
#   trua1b_2_<temp>.dependent-sites.tsv   same for allele trua1b_2
#   trua1b-<temp>.dependent.tsv           final dependent sites after the two-allele AND (intersection)
#
# Plot-table construction criteria (as in the paper's dependent_site_loss_effect):
#   - background = all WT-positive sites at that temperature (the WT_<temp> positive
#     set of the upstream Fisher step);
#   - WT_Ratio from the per-allele tables; Mut_Ratio = mean of the two alleles'
#     Mut_Ratio (verified to match, site by site, the Mut_Ratio carried in
#     dependent.tsv);
#   - site_status = dependent (present in the final intersection table) / wt_only
#     (the rest).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly TEMPERATURE_BB="$ROOT_DIR/Temperature/STAT/dependent_betabinomial"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# Join the three tables into one plot table. Row order follows line1; line1 and line2
# are verified to cover the identical site set (differing only in each allele's
# Mut_Ratio), and dependent is a subset of it.
write_plot_table() {
  local line1_tsv="$1"
  local line2_tsv="$2"
  local dependent_tsv="$3"
  local temp="$4"
  local output_tsv="$5"
  local rna_label="$6"

  local n_line1 n_line2
  n_line1=$(( $(wc -l < "$line1_tsv") - 1 ))
  n_line2=$(( $(wc -l < "$line2_tsv") - 1 ))
  [[ "$n_line1" -eq "$n_line2" ]] || {
    echo "[data_prep] ${rna_label}/${temp}: line1 (${n_line1}) and line2 (${n_line2}) site counts differ" >&2
    exit 1
  }

  awk -F'\t' -v OFS='\t' -v temp="$temp" -v rna_label="$rna_label" '
    BEGIN { print "Temperature", "chr", "pos", "strand", "context", "site_id",
                  "WT_Ratio", "Mut_Ratio", "Delta_Ratio", "site_status" }
    ARGIND == 1 {
      if (FNR == 1) { next }
      key = $1 SUBSEP $2 SUBSEP $3
      if (!(key in seen)) { seen[key] = 1; order[++n] = key }
      chr[key] = $1; pos[key] = $2; strand[key] = $3; context[key] = $4
      wt_ratio[key] = $5; mut_ratio_1[key] = $7
      next
    }
    ARGIND == 2 {
      if (FNR == 1) { next }
      key = $1 SUBSEP $2 SUBSEP $3
      mut_ratio_2[key] = $7
      next
    }
    ARGIND == 3 {
      if (FNR == 1) { next }
      dependent[$1 SUBSEP $2 SUBSEP $3] = 1
      next
    }
    END {
      for (i = 1; i <= n; i++) {
        key = order[i]
        if (!(key in mut_ratio_2)) { n_missing++ ; continue }
        mut_ratio = (mut_ratio_1[key] + mut_ratio_2[key]) / 2
        delta = wt_ratio[key] - mut_ratio
        status = (key in dependent) ? "dependent" : "wt_only"
        if (status == "dependent") { n_dependent++ }
        print temp, chr[key], pos[key], strand[key], context[key],
              (chr[key] "|" pos[key] "|" strand[key]),
              wt_ratio[key], mut_ratio, delta, status
      }
      printf "[data_prep] %s/%s: %d WT sites, %d dependent%s\n",
             rna_label, temp, n, n_dependent + 0,
             (n_missing + 0 > 0 ? sprintf(" (%d sites missing allele records, skipped)", n_missing) : "") > "/dev/stderr"
    }
  ' "$line1_tsv" "$line2_tsv" "$dependent_tsv" > "$output_tsv"
}

# Directory prefix of one RNA type × temperature combination
write_combo() {
  local rna_type="$1"   # mRNA / rRNA
  local temp="$2"
  local src_dir="${TEMPERATURE_BB}/${rna_type}/${temp}/trua1b"

  local line1_tsv="${src_dir}/trua1b_1_${temp}.dependent-sites.tsv"
  local line2_tsv="${src_dir}/trua1b_2_${temp}.dependent-sites.tsv"
  local dependent_tsv="${src_dir}/trua1b-${temp}.dependent.tsv"
  local output_tsv="${SCRIPT_DIR}/${rna_type}_${temp}_dependent_site_loss_effect.tsv"

  require_file "$line1_tsv"
  require_file "$line2_tsv"
  require_file "$dependent_tsv"

  write_plot_table "$line1_tsv" "$line2_tsv" "$dependent_tsv" "$temp" \
    "$output_tsv" "$rna_type"
}

main() {
  local rna_type temp
  for rna_type in mRNA rRNA; do
    for temp in 21 27; do
      write_combo "$rna_type" "$temp"
    done
  done
}

main "$@"
