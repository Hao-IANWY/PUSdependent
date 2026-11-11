#!/usr/bin/env bash
# data_prep.sh — Aggregate adjustGap modification fractions over all sites: WT + 7 mutants × 3 tissues
# Purpose: build the plot.R data interface merged_sample_ratios.tsv
# (tissue, sample, chr, pos, strand, ratio).
#
# Criteria (replicate set strictly aligned with the one used for dependent_betabinomial
# site calling):
#   WT     = mean(ratio) over the 3 replicates Treat-WT-<tissue>.rep{1,2,3}.adjusted.tsv.
#            rep4 is a cross-batch count-only replicate, excluded from site calling and
#            likewise excluded here; spot checks confirm this mean equals the
#            Treat_avg_ratio of fisher WT_<Tissue>_sites.Fisher.tsv.
#   Mutant = mean over the 3 replicates Treat-<g>_<line>-<tissue>.rep{1,2,3} per allele,
#            then averaged over the two alleles — the same computation chain as the
#            <g>_1/<g>_2/Mut_Ratio columns of the dependent tables.
#   WT ratios are computed from adjusted.tsv rather than Fisher.tsv because the tRNA
#   Fisher.tsv is a filtered subset (Root keeps only 623 sites) and cannot serve as a
#   full-site backfill source.
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
readonly ADJUST_DIR="${ROOT_DIR}/tRNA/STAT/adjustGap/ncRNA"
readonly GENOTYPES=(trua1a trua1b trua2 trua3 trua4 trua5 trua6)
readonly TISSUES=(Root leaf Flower)
readonly OUTPUT_TSV="${SCRIPT_DIR}/merged_sample_ratios.tsv"
readonly TMP_DIR="${SCRIPT_DIR}/.merge_heatmap_tmp"

cleanup() {
  rm -rf "${TMP_DIR}"
}
trap cleanup EXIT

require_nonempty_file() {
  local file="$1"
  [[ -s "${file}" ]] || {
    echo "[data_prep] Missing or empty input file: ${file}" >&2
    exit 1
  }
}

# adjusted.tsv columns: chr pos strand context ratio depth gap (ratio is always numeric; ratio=0 when depth=0)
mean_ratio() {
  # $1=output, remaining args=adjusted.tsv files of one sample; outputs chr pos strand mean(ratio)
  local out="$1"
  shift
  local f
  for f in "$@"; do
    require_nonempty_file "${f}"
  done
  awk -F '\t' -v OFS='\t' '
    FNR == 1 { next }
    { key = $1 SUBSEP $2 SUBSEP $3; sum[key] += $5; n[key]++ }
    END {
      for (k in sum) {
        split(k, a, SUBSEP)
        printf "%s\t%s\t%s\t%.10g\n", a[1], a[2], a[3], sum[k] / n[k]
      }
    }
  ' "$@" > "${out}"
}

merge_two_lines() {
  # $1=allele-1 means $2=allele-2 means $3=output; (v1+v2)/2, keep the single value when one side is missing
  awk -F '\t' -v OFS='\t' '
    NR == FNR { v[$1 SUBSEP $2 SUBSEP $3] = $4; next }
    {
      k = $1 SUBSEP $2 SUBSEP $3
      if (k in v) { printf "%s\t%s\t%s\t%.10g\n", $1, $2, $3, (v[k] + $4) / 2; delete v[k] }
      else        { printf "%s\t%s\t%s\t%.10g\n", $1, $2, $3, $4 }
    }
    END {
      for (k in v) {
        split(k, a, SUBSEP)
        printf "%s\t%s\t%s\t%.10g\n", a[1], a[2], a[3], v[k]
      }
    }
  ' "$1" "$2" > "$3"
}

append_with_keys() {
  # $1=tissue $2=sample $3=mean-ratio file -> appended to OUTPUT_TSV
  awk -F '\t' -v OFS='\t' -v tissue="$1" -v sample="$2" \
    '{ print tissue, sample, $0 }' "$3" >> "${OUTPUT_TSV}"
}

main() {
  mkdir -p "${TMP_DIR}"
  printf 'tissue\tsample\tchr\tpos\tstrand\tratio\n' > "${OUTPUT_TSV}"

  local tissue genotype
  for tissue in "${TISSUES[@]}"; do
    local wt_mean="${TMP_DIR}/WT.${tissue}.tsv"
    mean_ratio "${wt_mean}" \
      "${ADJUST_DIR}/Treat-WT-${tissue}.rep1.adjusted.tsv" \
      "${ADJUST_DIR}/Treat-WT-${tissue}.rep2.adjusted.tsv" \
      "${ADJUST_DIR}/Treat-WT-${tissue}.rep3.adjusted.tsv"
    append_with_keys "${tissue}" "WT" "${wt_mean}"
    echo "[data_prep] ${tissue} WT: $(wc -l < "${wt_mean}") sites"

    for genotype in "${GENOTYPES[@]}"; do
      local line1_mean="${TMP_DIR}/${genotype}_1.${tissue}.tsv"
      local line2_mean="${TMP_DIR}/${genotype}_2.${tissue}.tsv"
      local merged="${TMP_DIR}/${genotype}.${tissue}.tsv"
      mean_ratio "${line1_mean}" \
        "${ADJUST_DIR}/Treat-${genotype}_1-${tissue}.rep1.adjusted.tsv" \
        "${ADJUST_DIR}/Treat-${genotype}_1-${tissue}.rep2.adjusted.tsv" \
        "${ADJUST_DIR}/Treat-${genotype}_1-${tissue}.rep3.adjusted.tsv"
      mean_ratio "${line2_mean}" \
        "${ADJUST_DIR}/Treat-${genotype}_2-${tissue}.rep1.adjusted.tsv" \
        "${ADJUST_DIR}/Treat-${genotype}_2-${tissue}.rep2.adjusted.tsv" \
        "${ADJUST_DIR}/Treat-${genotype}_2-${tissue}.rep3.adjusted.tsv"
      merge_two_lines "${line1_mean}" "${line2_mean}" "${merged}"
      append_with_keys "${tissue}" "${genotype}" "${merged}"
      echo "[data_prep] ${tissue} ${genotype}: $(wc -l < "${merged}") sites"
    done
  done

  local row_count
  row_count="$(( $(wc -l < "${OUTPUT_TSV}") - 1 ))"
  echo "[data_prep] Generated: ${OUTPUT_TSV} (${row_count} rows)"
}

main "$@"
