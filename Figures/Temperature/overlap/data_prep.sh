#!/usr/bin/env bash
# data_prep.sh — 「27℃ 温度上调位点」与「27℃ trua1b dependent 位点」的 overlap 表
#               （mRNA / rRNA / tRNA，各一张）
#
# 职责：只产 TSV，不画图。
#
# 两套集合各三个 RNA 类型，来源：
#   27℃ 上调        : Temperature/STAT/Temperature_response/<type>/WT_27_vs_WT_21.response.tsv
#                     （direction == up_in_WT_27；tRNA 侧走 Temperature_tRNA 的 ncRNA）
#   27℃ dependent   : Temperature/STAT/dependent_methylKit/<type>/27/trua1b/trua1b-27.dependent.tsv
#
# 位点主键 = chr + pos + strand；两个集合取并集逐位点打标记，输出 combination
# （up_only / dependent_only / both），交给 plot.R 数三种计数画 2-set Venn。
#
# tRNA 侧注意：Temperature_tRNA 的参考是 tRNA-unique_sn_snoRNA.fa，混有 sn/snoRNA，
# 两个集合都必须按 ^tRNA- 过滤（与 Figures/Temperature 其余图组同口径）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"

readonly RESPONSE_ROOT="$ROOT_DIR/Temperature/STAT/Temperature_response"
readonly RESPONSE_TRNA="$ROOT_DIR/Temperature_tRNA/STAT/Temperature_response/ncRNA"
readonly DMK_ROOT="$ROOT_DIR/Temperature/STAT/dependent_methylKit"
readonly DMK_TRNA="$ROOT_DIR/Temperature_tRNA/STAT/dependent_methylKit"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

# trna_only 用 0/1 传递：awk 里非空字符串（"false"）恒为真，不能用 true/false
write_overlap_table() {
  local up_tsv="$1"
  local dep_tsv="$2"
  local trna_only="$3"
  local output_tsv="$4"
  local rna_label="$5"

  awk -F'\t' -v OFS='\t' -v trna_only="$trna_only" -v rna_label="$rna_label" '
    function keep() { return !(trna_only == 1 && $1 !~ /^tRNA-/) }
    function register(k, c, p, s, ctx) {
      if (!(k in seen)) { seen[k] = 1; order[++n] = k }
      chr[k] = c; pos[k] = p; strand[k] = s
      if (ctx != "" && !(k in context)) { context[k] = ctx }
    }
    BEGIN { print "chr", "pos", "strand", "context", "site_id", "up_in_WT_27",
                  "trua1b_dependent_27", "combination" }
    ARGIND == 1 {
      if (FNR == 1 || !keep() || $9 != "up_in_WT_27") { next }
      gsub(/\r/, "")
      k = $1 SUBSEP $2 SUBSEP $3
      register(k, $1, $2, $3, $4)
      up[k] = 1
      next
    }
    ARGIND == 2 {
      if (FNR == 1 || !keep()) { next }
      gsub(/\r/, "")
      k = $1 SUBSEP $2 SUBSEP $3
      register(k, $1, $2, $3, $4)
      dep[k] = 1
      next
    }
    END {
      for (i = 1; i <= n; i++) {
        k = order[i]
        u = (k in up) ? 1 : 0
        d = (k in dep) ? 1 : 0
        combo = (u && d) ? "both" : ((u) ? "up_only" : "dependent_only")
        count[combo]++
        print chr[k], pos[k], strand[k], context[k], (chr[k] "|" pos[k] "|" strand[k]),
              u, d, combo
      }
      printf "[data_prep] %s: up_only=%d, dependent_only=%d, both=%d\n",
             rna_label, count["up_only"] + 0, count["dependent_only"] + 0,
             count["both"] + 0 > "/dev/stderr"
    }
  ' "$up_tsv" "$dep_tsv" > "$output_tsv"
}

write_pair() {
  local rna_type="$1"
  local trna_only="$2"
  local up_tsv dep_tsv

  if [[ "$rna_type" == "tRNA" ]]; then
    up_tsv="${RESPONSE_TRNA}/WT_27_vs_WT_21.response.tsv"
    dep_tsv="${DMK_TRNA}/ncRNA/27/trua1b/trua1b-27.dependent.tsv"
  else
    up_tsv="${RESPONSE_ROOT}/${rna_type}/WT_27_vs_WT_21.response.tsv"
    dep_tsv="${DMK_ROOT}/${rna_type}/27/trua1b/trua1b-27.dependent.tsv"
  fi

  require_file "$up_tsv"
  require_file "$dep_tsv"

  write_overlap_table "$up_tsv" "$dep_tsv" "$trna_only" \
    "${SCRIPT_DIR}/${rna_type}_27up_vs_trua1b_dependent_overlap.tsv" "$rna_type"
}

main() {
  write_pair mRNA 0
  write_pair rRNA 0
  write_pair tRNA 1
}

main "$@"
