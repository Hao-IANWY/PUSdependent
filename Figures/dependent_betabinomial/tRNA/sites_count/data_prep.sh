#!/usr/bin/env bash
# data_prep.sh — Count writer-dependent Ψ sites per genotype × tissue from
# tRNA/STAT/dependent_betabinomial/ncRNA
# Purpose: build the TSV data file read by plot.R; no plotting here.
#
# Scope: count only tRNA-derived sites. The ncRNA reference of the tRNA workflow is
# tRNA-unique_sn_snoRNA.fa (589 sequences: 220 tRNAs "tRNA-aa-codon-n" + 287 snoRNAs
# "AT...-sno" + 82 snRNAs "AT...-sn"); the dependent tables also contain sn/snoRNA
# sites, which are dropped by the chr prefix filter.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
INPUT_DIR="$ROOT_DIR/tRNA/STAT/dependent_betabinomial/ncRNA"
OUT_TSV="$SCRIPT_DIR/dependent_site_count.tsv"

readonly TRNA_CHR_PATTERN='^tRNA-'

tissues=("Root" "leaf" "Flower")
genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")

printf 'genotype\ttissue\tcount\n' > "$OUT_TSV"

for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    dep_file="$INPUT_DIR/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
    if [[ -s "$dep_file" ]]; then
      # One dependent site per row, header first; count only rows whose chr starts with tRNA-
      count=$(awk -F'\t' -v pat="$TRNA_CHR_PATTERN" 'NR > 1 && $1 ~ pat {n++} END {print n+0}' "$dep_file")
    else
      count=0
    fi
    printf '%s\t%s\t%s\n' "$genotype" "$tissue" "$count" >> "$OUT_TSV"
  done
done

echo "[data_prep] Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) data rows)"
