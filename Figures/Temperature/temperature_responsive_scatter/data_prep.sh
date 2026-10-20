#!/usr/bin/env bash
# data_prep.sh — 汇总 21/27 度 WT 温度响应位点的散点图数据表（mRNA / rRNA / tRNA 各一张）
#
# 职责：只产 TSV，不画图。三张表同构，plot.R 逐张出图。
#
# 数据来源（两套 STAT 流程，产物列名完全一致）：
#   mRNA / rRNA : Temperature/STAT/Temperature_response/<type>/WT_27_vs_WT_21.all-sites.tsv
#   tRNA        : Temperature_tRNA/STAT/Temperature_response/ncRNA/WT_27_vs_WT_21.all-sites.tsv
#
# all-sites.tsv 就是候选位点全表（= 两个温度 fisher 阳性位点并集，逐位点重算
# 两温度修饰比与 Mock 平均深度），与 archive 论文版 temperature_responsive_sites
# 里由 data.sh 现拼的 temperature_responsive_scatter.tsv 语义等价，故这里不再
# 回原始 fisher/adjustGap 重算，直接取库内产物，列名做一次规范化。
#
# tRNA 侧注意：Temperature_tRNA 的参考是 tRNA-unique_sn_snoRNA.fa，混有 sn/snoRNA，
# 必须按 ^tRNA- 过滤，只保留 tRNA 部分（与 Figures/dependent_methylKit/tRNA 同口径）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

readonly TEMPERATURE_ALL_SITES="$ROOT_DIR/Temperature/STAT/Temperature_response"
readonly TEMPERATURE_TRNA_ALL_SITES="$ROOT_DIR/Temperature_tRNA/STAT/Temperature_response/ncRNA/WT_27_vs_WT_21.all-sites.tsv"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# 取 all-sites.tsv -> 画图表。
# 保留列：chr pos strand context WT_21_Ratio WT_27_Ratio Delta_Ratio_27_minus_21
#         Mock_Mean_Depth_WT_21 Mock_Mean_Depth_WT_27 p_value q_value regulation
# regulation 由上游 direction 列改写而来：空值 -> non-responsive，
# up_in_WT_27 -> up-regulated（27 度修饰升高），down_in_WT_27 -> down-regulated。
# 判定口径见各 STAT/Temperature_response/run.sh：methylKit F 检验、
# SLIM q < 0.05（双侧 FDR）、|Delta| >= 0.05、两个温度组每个重复 depth >= 10。
write_plot_table() {
  local input_tsv="$1"
  local output_tsv="$2"
  local trna_only="$3"

  # trna_only 用 0/1 传递：awk 里非空字符串（"false"）恒为真，不能用 true/false
  awk -F'\t' -v OFS='\t' -v trna_only="$trna_only" '
    BEGIN { print "chr", "pos", "strand", "context", "WT_21_Ratio", "WT_27_Ratio",
                  "Delta_Ratio_27_minus_21", "Mock_Mean_Depth_WT_21", "Mock_Mean_Depth_WT_27",
                  "p_value", "q_value", "regulation" }
    NR == 1 { next }
    {
      if (trna_only == 1 && $1 !~ /^tRNA-/) { next }
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
    write_plot_table "$input_tsv" "$output_tsv" "0"
  done

  # tRNA：从 ncRNA 全表里按 ^tRNA- 过滤
  require_file "$TEMPERATURE_TRNA_ALL_SITES"
  write_plot_table "$TEMPERATURE_TRNA_ALL_SITES" \
    "${SCRIPT_DIR}/tRNA_temperature_responsive_scatter.tsv" "1"
}

main "$@"
