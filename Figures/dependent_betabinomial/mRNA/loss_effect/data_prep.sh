#!/usr/bin/env bash
# data_prep.sh — Assemble the loss-of-modification data for mRNA writer-dependent Ψ
# sites across 7 trua mutants × 3 tissues.
# Purpose: build dependent_site_loss_effect_dotplot.tsv for plot.R; no plotting here.
#
# One row per WT-positive mRNA site (the full background set), labeled by whether it
# belongs to the dependent set of the corresponding genotype × tissue:
#   WT_Ratio  from dependent-sites.tsv WT_Ratio (mean Treat modification fraction of the
#             WT-positive candidates)
#   Mut_Ratio from dependent-sites.tsv Mut_Ratio (mean of the unfiltered adjustGap pileup
#             over the 3 Treat replicates), then averaged over the two alleles
#
# Note: Mut_Ratio must come from the unfiltered adjustGap criterion, not from the Fisher
# tables. Fisher filtering zeroes weak residual mutant signals that are still present
# (e.g. AT1G04210.1 has adjustGap ratios 0.000/0.009174/0.012658 -> mean 0.007277, yet
# the Fisher table drops it and reports Mut_Ratio=0). The Mut_Ratio in
# dependent-sites.tsv is exactly the mean over the 3 Treat replicates computed from
# adjustGap by call_dependent.py.
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
DEPENDENT_DIR="$ROOT_DIR/polyARNA/STAT/dependent_betabinomial/mRNA"
OUT_TSV="$SCRIPT_DIR/dependent_site_loss_effect_dotplot.tsv"

genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")
tissues=("Root" "leaf" "Flower")

printf 'panel_id\tgenotype\ttissue\tchr\tpos\tstrand\tcontext\tsite_id\tWT_Ratio\tMut_Ratio\tDelta_Ratio\tsite_status\n' > "$OUT_TSV"

for genotype in "${genotypes[@]}"; do
  for tissue in "${tissues[@]}"; do
    line1_file="$DEPENDENT_DIR/$tissue/$genotype/${genotype}_1-${tissue}.dependent-sites.tsv"
    line2_file="$DEPENDENT_DIR/$tissue/$genotype/${genotype}_2-${tissue}.dependent-sites.tsv"

    for f in "$line1_file" "$line2_file"; do
      [[ -s "$f" ]] || { echo "[data_prep] Missing input file: $f" >&2; exit 1; }
    done

    awk -F '\t' -v OFS='\t' \
        -v genotype="$genotype" \
        -v tissue="$tissue" \
        -v panel_id="${genotype}_${tissue}" '
        function validate_required(required_names,    i, name) {
            split(required_names, fields, ",")
            for (i in fields) {
                name = fields[i]
                if (!(name in col)) {
                    printf("Required column missing in %s: %s\n", FILENAME, name) > "/dev/stderr"
                    exit 1
                }
            }
        }
        FNR == 1 {
            file_index += 1
            delete col
            for (i = 1; i <= NF; i++) {
                sub(/\r$/, "", $i)
                col[$i] = i
            }
            validate_required("chr,pos,strand,context,WT_Ratio,Mut_Ratio,is_dependent")
            next
        }
        {
            chr = $(col["chr"])
            pos = $(col["pos"])
            strand = $(col["strand"])
            site_id = chr "|" pos "|" strand
            wt_ratio = $(col["WT_Ratio"])
            mut_ratio = $(col["Mut_Ratio"])
            dep_raw = $(col["is_dependent"])
            sub(/\r$/, "", dep_raw)
            is_dep = (tolower(dep_raw) == "true" ? 1 : 0)

            if (file_index == 1) {
                n += 1
                order[n] = site_id
                chr_map[site_id] = chr
                pos_map[site_id] = pos
                strand_map[site_id] = strand
                context_map[site_id] = $(col["context"])
                wt_map[site_id] = wt_ratio
                mut_sum[site_id] = mut_ratio
                mut_count[site_id] = 1
                dep_map[site_id] = is_dep
                next
            }

            # file_index == 2: merge with line1; accumulate Mut_Ratio; a site counts as
            # dependent only when is_dependent holds in BOTH alleles
            if (site_id in wt_map) {
                mut_sum[site_id] += mut_ratio
                mut_count[site_id] += 1
                dep_map[site_id] = (dep_map[site_id] && is_dep)
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                site_id = order[i]
                mut_ratio = (mut_count[site_id] > 0 ? mut_sum[site_id] / mut_count[site_id] : 0)
                delta = wt_map[site_id] - mut_ratio
                status = (dep_map[site_id] ? "dependent" : "wt_only")
                print panel_id, genotype, tissue, chr_map[site_id], pos_map[site_id], strand_map[site_id], \
                      context_map[site_id], site_id, wt_map[site_id], mut_ratio, delta, status
            }
        }
    ' "$line1_file" "$line2_file" >> "$OUT_TSV"

    echo "[data_prep] $genotype / $tissue done"
  done
done

echo "[data_prep] Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) data rows)"
