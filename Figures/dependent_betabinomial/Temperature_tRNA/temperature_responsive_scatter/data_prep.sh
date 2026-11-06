#!/usr/bin/env bash
# data_prep.sh — 汇总 21/27 度 WT 温度响应位点的散点图数据表（tRNA 一张）
#
# 职责：只产 TSV，不画图。mRNA / rRNA 版在姊妹目录
#   Figures/dependent_betabinomial/Temperature_mRNA/temperature_responsive_scatter/
#
# 数据来源：
#   tRNA : Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.all-sites.tsv
#
# all-sites.tsv 就是**候选位点全表**（两个温度 fisher 阳性位点并集，逐位点回
# adjustGap 重算两个温度的修饰比与 Mock 平均深度），与 archive 论文版由 data.sh
# 现拼的 temperature_responsive_scatter.tsv 语义等价，故这里不再回原始
# fisher/adjustGap 重算，直接取库内产物，列名做一次规范化。
#
# tRNA 侧注意：Temperature_tRNA 的参考是 tRNA-unique_sn_snoRNA.fa，混有 sn/snoRNA，
# 必须按 ^tRNA- 过滤，只保留 tRNA 部分（与 Figures/dependent_betabinomial/tRNA 同口径）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly TEMPERATURE_TRNA_ALL_SITES="$ROOT_DIR/Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.all-sites.tsv"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# 取 all-sites.tsv -> 画图表。
# 保留列：chr pos strand context WT_21_Ratio WT_27_Ratio Delta_Ratio_27_minus_21
#         Mock_Mean_Depth_WT_21 Mock_Mean_Depth_WT_27 p_value q_value regulation
# regulation 由上游 direction 列改写而来：空值 -> non-responsive，
# up_in_WT_27 -> up-regulated（27 度修饰升高），down_in_WT_27 -> down-regulated。
# 判定口径见 Temperature_tRNA/STAT/Temperature_response_betabinomial/run.sh：
# beta-binomial GLMM 全局 theta + 逐位点 LRT（--site-test lrt）、BH q < 0.05
# （双侧 FDR）、|Delta| >= 0.05、两个温度组每个重复 depth >= 10。
write_plot_table() {
  local input_tsv="$1"
  local output_tsv="$2"

  awk -F'\t' -v OFS='\t' '
    BEGIN { print "chr", "pos", "strand", "context", "WT_21_Ratio", "WT_27_Ratio",
                  "Delta_Ratio_27_minus_21", "Mock_Mean_Depth_WT_21", "Mock_Mean_Depth_WT_27",
                  "p_value", "q_value", "regulation" }
    NR == 1 { next }
    {
      if ($1 !~ /^tRNA-/) { next }
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
  require_file "$TEMPERATURE_TRNA_ALL_SITES"
  write_plot_table "$TEMPERATURE_TRNA_ALL_SITES" \
    "${SCRIPT_DIR}/tRNA_temperature_responsive_scatter.tsv"
}

main "$@"
