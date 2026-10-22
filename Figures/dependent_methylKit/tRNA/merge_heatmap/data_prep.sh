#!/usr/bin/env bash
# data_prep.sh — 聚合 adjustGap 全位点修饰比例：WT + 7 个突变体 × 3 组织
# 职责：生成 plot.R 的数据接口 merged_sample_ratios.tsv（tissue, sample, chr, pos, strand, ratio）。
#
# 口径（与 dependent_methylKit 位点鉴定所用重复集严格对齐）：
#   WT     = Treat-WT-<tissue>.rep{1,2,3}.adjusted.tsv 三重复 mean(ratio)。
#            rep4 是跨批次 count-only 重复、不进鉴定口径，这里同样排除；
#            已抽查验证该均值 == fisher WT_<Tissue>_sites.Fisher.tsv 的 Treat_avg_ratio。
#   突变体 = 每株系 Treat-<g>_<line>-<tissue>.rep{1,2,3} 三重复 mean → 两株系再 mean，
#            与 dependent 表 <g>_1/<g>_2/Mut_Ratio 的计算链一致。
#   archive 版 WT 比例取自 Fisher.tsv，但 tRNA 流程的 Fisher.tsv 是过滤子集（Root 仅 623 位点），
#   不能作全位点回填来源，故 WT 同样走 adjusted.tsv（与 archive 的唯一上游差异）。
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

# adjusted.tsv 列：chr pos strand context ratio depth gap（ratio 恒为数值，depth=0 时 ratio=0）
mean_ratio() {
  # $1=输出，其余=同一样本的若干 adjusted.tsv；输出 chr pos strand mean(ratio)
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
  # $1=株系1均值 $2=株系2均值 $3=输出；(v1+v2)/2，单侧缺失取单侧值
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
  # $1=tissue $2=sample $3=均值文件 → 追加到 OUTPUT_TSV
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
