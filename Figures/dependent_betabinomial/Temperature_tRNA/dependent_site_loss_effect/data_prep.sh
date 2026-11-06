#!/usr/bin/env bash
# data_prep.sh — 汇总 trua1b dependent 位点的修饰下降效应散点图数据表
#               （tRNA × 21℃ / 27℃，共 2 张）
#
# 职责：只产 TSV，不画图。mRNA / rRNA 版在姊妹目录
#   Figures/dependent_betabinomial/Temperature_mRNA/dependent_site_loss_effect/
#
# 数据来源：
#   tRNA : Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/<temp>/trua1b/
#
# 每个组合三张表：
#   trua1b_1_<temp>.dependent-sites.tsv   株系 trua1b_1 的逐位点表（含全部 WT 候选位点）
#   trua1b_2_<temp>.dependent-sites.tsv   株系 trua1b_2（原 trua1b_3，2026-11-04 正名）同上
#   trua1b-<temp>.dependent.tsv           两株系 AND（交集）后的最终 dependent 位点
#
# 画图表的构建口径（与 archive 论文版 dependent_site_loss_effect 一致）：
#   - 背景 = 该温度下的全部 WT 阳性位点（= 上游 fisher 的 WT_<temp> 阳性集）；
#   - WT_Ratio 取 per-line 表；Mut_Ratio = 两株系 Mut_Ratio 的均值
#     （已实测与 dependent.tsv 自带的 Mut_Ratio 逐位点完全一致）；
#   - site_status = dependent（落在最终交集表里的）/ wt_only（其余）。
#
# tRNA 侧注意：Temperature_tRNA 的参考是 tRNA-unique_sn_snoRNA.fa，混有 sn/snoRNA，
# 必须按 ^tRNA- 过滤，只保留 tRNA 部分（与 Figures/dependent_betabinomial/tRNA 同口径）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly TEMPERATURE_TRNA_BB="$ROOT_DIR/Temperature_tRNA/STAT/dependent_betabinomial"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# 三个表关联成一张画图表。order 取自 line1 的行序；已确认 line1 / line2 的位点集合
# 完全相同（只差各自株系的 Mut_Ratio），dependent 是其子集。
# 株系命名：trua1b_2 = 2026-11-04 起的新名（原 trua1b_3 正名，诺禾旧批 35 系同系）。
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
    echo "[data_prep] ${rna_label}/${temp}: line1 (${n_line1}) 与 line2 (${n_line2}) 位点数不一致" >&2
    exit 1
  }

  awk -F'\t' -v OFS='\t' -v temp="$temp" -v rna_label="$rna_label" '
    function keep() { return $1 ~ /^tRNA-/ }
    BEGIN { print "Temperature", "chr", "pos", "strand", "context", "site_id",
                  "WT_Ratio", "Mut_Ratio", "Delta_Ratio", "site_status" }
    ARGIND == 1 {
      if (FNR == 1 || !keep()) { next }
      key = $1 SUBSEP $2 SUBSEP $3
      if (!(key in seen)) { seen[key] = 1; order[++n] = key }
      chr[key] = $1; pos[key] = $2; strand[key] = $3; context[key] = $4
      wt_ratio[key] = $5; mut_ratio_1[key] = $7
      next
    }
    ARGIND == 2 {
      if (FNR == 1 || !keep()) { next }
      key = $1 SUBSEP $2 SUBSEP $3
      mut_ratio_2[key] = $7
      next
    }
    ARGIND == 3 {
      if (FNR == 1 || !keep()) { next }
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
             (n_missing + 0 > 0 ? sprintf(" (%d 位点缺株系记录,已跳过)", n_missing) : "") > "/dev/stderr"
    }
  ' "$line1_tsv" "$line2_tsv" "$dependent_tsv" > "$output_tsv"
}

# 单个温度的组合
write_combo() {
  local temp="$1"
  local src_dir="${TEMPERATURE_TRNA_BB}/ncRNA/${temp}/trua1b"

  local line1_tsv="${src_dir}/trua1b_1_${temp}.dependent-sites.tsv"
  local line2_tsv="${src_dir}/trua1b_2_${temp}.dependent-sites.tsv"
  local dependent_tsv="${src_dir}/trua1b-${temp}.dependent.tsv"
  local output_tsv="${SCRIPT_DIR}/tRNA_${temp}_dependent_site_loss_effect.tsv"

  require_file "$line1_tsv"
  require_file "$line2_tsv"
  require_file "$dependent_tsv"

  write_plot_table "$line1_tsv" "$line2_tsv" "$dependent_tsv" "$temp" \
    "$output_tsv" "tRNA"
}

main() {
  local temp
  for temp in 21 27; do
    write_combo "$temp"
  done
}

main "$@"
