#!/usr/bin/env bash
# data_prep.sh — 汇总 trua1b dependent 位点的修饰下降效应散点图数据表
#               （mRNA / rRNA / tRNA × 21℃ / 27℃，共 6 张）
#
# 职责：只产 TSV，不画图。
#
# 数据来源（Temperature_tRNA 与 Temperature 同构，仅 mpileup 参考不同）：
#   mRNA / rRNA : Temperature/STAT/dependent_betabinomial/<type>/<temp>/trua1b/
#   tRNA        : Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/<temp>/trua1b/
#
# 每个组合三张表：
#   trua1b_1_<temp>.dependent-sites.tsv   株系 trua1b_1 的逐位点表（含全部 WT 候选位点）
#   trua1b_3_<temp>.dependent-sites.tsv   株系 trua1b_3 同上
#   trua1b-<temp>.dependent.tsv           两株系 AND（交集）后的最终 dependent 位点
#
# 画图表的构建口径（与 archive 论文版 dependent_site_loss_effect 一致）：
#   - 背景 = 该温度下的全部 WT 阳性位点（= 上游 fisher 的 WT_<temp> 阳性集）；
#   - WT_Ratio 取 per-line 表；Mut_Ratio = 两株系 Mut_Ratio 的均值
#     （已实测与 dependent.tsv 自带的 Mut_Ratio 逐位点完全一致）；
#   - site_status = dependent（落在最终交集表里的）/ wt_only（其余）。
#
# tRNA 侧注意：Temperature_tRNA 的参考是 tRNA-unique_sn_snoRNA.fa，混有 sn/snoRNA，
# 必须按 ^tRNA- 过滤，只保留 tRNA 部分（与 Figures/dependent_methylKit/tRNA 同口径）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly TEMPERATURE_BB="$ROOT_DIR/Temperature/STAT/dependent_betabinomial"
readonly TEMPERATURE_TRNA_BB="$ROOT_DIR/Temperature_tRNA/STAT/dependent_betabinomial"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# 三个表关联成一张画图表。order 取自 line1 的行序；已确认 line1 / line3 的位点集合
# 完全相同（只差各自株系的 Mut_Ratio），dependent 是其子集。
# trna_only 用 0/1 传递：awk 里非空字符串（"false"）恒为真，不能用 true/false。
write_plot_table() {
  local line1_tsv="$1"
  local line3_tsv="$2"
  local dependent_tsv="$3"
  local temp="$4"
  local trna_only="$5"
  local output_tsv="$6"
  local rna_label="$7"

  local n_line1 n_line3
  n_line1=$(( $(wc -l < "$line1_tsv") - 1 ))
  n_line3=$(( $(wc -l < "$line3_tsv") - 1 ))
  [[ "$n_line1" -eq "$n_line3" ]] || {
    echo "[data_prep] ${rna_label}/${temp}: line1 (${n_line1}) 与 line3 (${n_line3}) 位点数不一致" >&2
    exit 1
  }

  awk -F'\t' -v OFS='\t' -v temp="$temp" -v trna_only="$trna_only" -v rna_label="$rna_label" '
    function keep() { return !(trna_only == 1 && $1 !~ /^tRNA-/) }
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
      mut_ratio_3[key] = $7
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
        if (!(key in mut_ratio_3)) { n_missing++ ; continue }
        mut_ratio = (mut_ratio_1[key] + mut_ratio_3[key]) / 2
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
  ' "$line1_tsv" "$line3_tsv" "$dependent_tsv" > "$output_tsv"
}

# 单个 RNA 类型 × 温度的目录前缀
write_combo() {
  local rna_type="$1"   # mRNA / rRNA / tRNA
  local temp="$2"
  local trna_only="$3"
  local src_dir

  if [[ "$rna_type" == "tRNA" ]]; then
    src_dir="${TEMPERATURE_TRNA_BB}/ncRNA/${temp}/trua1b"
  else
    src_dir="${TEMPERATURE_BB}/${rna_type}/${temp}/trua1b"
  fi

  local line1_tsv="${src_dir}/trua1b_1_${temp}.dependent-sites.tsv"
  local line3_tsv="${src_dir}/trua1b_3_${temp}.dependent-sites.tsv"
  local dependent_tsv="${src_dir}/trua1b-${temp}.dependent.tsv"
  local output_tsv="${SCRIPT_DIR}/${rna_type}_${temp}_dependent_site_loss_effect.tsv"

  require_file "$line1_tsv"
  require_file "$line3_tsv"
  require_file "$dependent_tsv"

  write_plot_table "$line1_tsv" "$line3_tsv" "$dependent_tsv" "$temp" "$trna_only" \
    "$output_tsv" "$rna_type"
}

main() {
  local rna_type temp
  for rna_type in mRNA rRNA tRNA; do
    for temp in 21 27; do
      if [[ "$rna_type" == "tRNA" ]]; then
        write_combo "$rna_type" "$temp" 1
      else
        write_combo "$rna_type" "$temp" 0
      fi
    done
  done
}

main "$@"
