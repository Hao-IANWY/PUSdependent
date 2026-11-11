#!/usr/bin/env bash
# data_prep.sh — Build the WT-ratio background and mark TRUA-dependent sites.
# Inputs:  polyARNA/STAT/dependent_betabinomial/mRNA/{Root,leaf,Flower}/<genotype>/
# Outputs: dependent_site_ratio_distribution_sites.tsv
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
INPUT_DIR="$ROOT_DIR/polyARNA/STAT/dependent_betabinomial/mRNA"
OUT_TSV="$SCRIPT_DIR/dependent_site_ratio_distribution_sites.tsv"

tissues=("Root" "leaf" "Flower")
genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT

dependent_map="$tmpdir/dependent_map.tsv"
dependent_union="$tmpdir/dependent_union.tsv"
printf 'Tissue\tchr\tpos\tstrand\tGenotype\n' > "$dependent_map"
printf 'Tissue\tchr\tpos\tstrand\tDependent_Genotypes\n' > "$dependent_union"

# Collect genotype-level dependent calls, then collapse duplicate sites to a
# comma-separated genotype list.
for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    dep_file="$INPUT_DIR/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
    if [[ ! -s "$dep_file" ]]; then
      printf '[data_prep] Missing dependent file: %s\n' "$dep_file" >&2
      exit 1
    fi
    awk -F '\t' -v OFS='\t' -v tissue="$tissue" -v genotype="$genotype" '
      NR == 1 { next }
      NF >= 3 { print tissue, $1, $2, $3, genotype }
    ' "$dep_file" >> "$dependent_map"
  done
done

awk -F '\t' -v OFS='\t' '
  NR == 1 { next }
  {
    key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
    if (!(key in seen)) {
      seen[key] = 1
      order[++n] = key
      genotypes[key] = $5
    } else if (("," genotypes[key] ",") !~ ("," $5 ",")) {
      genotypes[key] = genotypes[key] "," $5
    }
  }
  END {
    for (i = 1; i <= n; i++) {
      key = order[i]
      split(key, field, SUBSEP)
      print field[1], field[2], field[3], field[4], genotypes[key]
    }
  }
' "$dependent_map" >> "$dependent_union"

printf '[data_prep] Dependent union: %d sites\n' "$(( $(wc -l < "$dependent_union") - 1 ))"
printf 'Tissue\tchr\tpos\tstrand\tcontext\tWT_Ratio\tSite_ID\tSite_Status\tDependent_Genotypes\n' > "$OUT_TSV.tmp"

# All genotype-specific WT tables are copies of the same tissue-level WT set.
# Use trua1a as the canonical copy and verify that the other copies agree.
for tissue in "${tissues[@]}"; do
  wt_file="$INPUT_DIR/$tissue/trua1a/WT-${tissue}_sites.Treat_avg_ratio.tsv"
  if [[ ! -s "$wt_file" ]]; then
    printf '[data_prep] Missing WT table: %s\n' "$wt_file" >&2
    exit 1
  fi
  for genotype in "${genotypes[@]:1}"; do
    other_wt_file="$INPUT_DIR/$tissue/$genotype/WT-${tissue}_sites.Treat_avg_ratio.tsv"
    if [[ ! -s "$other_wt_file" ]]; then
      printf '[data_prep] Missing WT table: %s\n' "$other_wt_file" >&2
      exit 1
    fi
    if ! cmp -s "$wt_file" "$other_wt_file"; then
      printf '[data_prep] WT tables differ for %s: %s vs %s\n' "$tissue" "$wt_file" "$other_wt_file" >&2
      exit 1
    fi
  done

  awk -F '\t' -v OFS='\t' -v tissue="$tissue" -v dep_file="$dependent_union" '
    FILENAME == dep_file {
      if (FNR == 1) next
      if ($1 != tissue) next
      key = $1 SUBSEP $2 SUBSEP $3 SUBSEP $4
      dep[key] = $5
      expected_dependent++
      next
    }
    FNR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      required = "chr pos strand context Treat_avg_ratio Pass"
      n = split(required, req, " ")
      for (i = 1; i <= n; i++) {
        if (!(req[i] in col)) {
          printf("Required column missing in %s: %s\n", FILENAME, req[i]) > "/dev/stderr"
          exit 1
        }
      }
      next
    }
    {
      if ($(col["Pass"]) != 1) next
      wt_ratio = $(col["Treat_avg_ratio"]) + 0
      if (wt_ratio < 0 || wt_ratio > 1) next

      chr = $(col["chr"])
      pos = $(col["pos"])
      strand = $(col["strand"])
      context = $(col["context"])
      key = tissue SUBSEP chr SUBSEP pos SUBSEP strand
      site_id = chr ":" pos ":" strand
      status = (key in dep ? "dependent" : "wt_only")
      dep_genotypes = (key in dep ? dep[key] : "NA")
      if (key in dep) matched_dependent++
      print tissue, chr, pos, strand, context, wt_ratio, site_id, status, dep_genotypes
    }
    END {
      if (matched_dependent != expected_dependent) {
        printf("Dependent/WT site mismatch for %s: expected %d, matched %d\n", tissue, expected_dependent, matched_dependent) > "/dev/stderr"
        exit 1
      }
    }
  ' "$dependent_union" "$wt_file" >> "$OUT_TSV.tmp"
done

mv "$OUT_TSV.tmp" "$OUT_TSV"
printf '[data_prep] Generated: %s (%d data rows)\n' "$OUT_TSV" "$(( $(wc -l < "$OUT_TSV") - 1 ))"
