#!/usr/bin/env bash
# data_prep.sh — 「27℃ 温度上调位点」与「27℃ trua1b dependent 位点」的 overlap 表
#               （mRNA / rRNA，各一张）
#
# 职责：只产 TSV，不画图。tRNA 版在姊妹目录
#   Figures/dependent_betabinomial/Temperature_tRNA/overlap/
#
# 两个集合各两个 RNA 类型，来源：
#   27℃ 上调        : Temperature/STAT/Temperature_response_betabinomial/<type>/WT_27_vs_WT_21.response.tsv
#                      （direction == up_in_WT_27）
#   27℃ dependent   : Temperature/STAT/dependent_betabinomial/<type>/27/trua1b/trua1b-27.dependent.tsv
#
# 位点主键 = chr + pos + strand；两个集合取并集逐位点打标记，输出 combination
# （up_only / dependent_only / both），交给 plot.R 数三种计数画 2-set Venn。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"

readonly RESPONSE_ROOT="$ROOT_DIR/Temperature/STAT/Temperature_response_betabinomial"
readonly BB_ROOT="$ROOT_DIR/Temperature/STAT/dependent_betabinomial"

require_file() {
  [[ -s "$1" ]] || { echo "[data_prep] Missing or empty input: $1" >&2; exit 1; }
}

write_overlap_table() {
  local up_tsv="$1"
  local dep_tsv="$2"
  local output_tsv="$3"
  local rna_label="$4"

  awk -F'\t' -v OFS='\t' -v rna_label="$rna_label" '
    function register(k, c, p, s, ctx) {
      if (!(k in seen)) { seen[k] = 1; order[++n] = k }
      chr[k] = c; pos[k] = p; strand[k] = s
      if (ctx != "" && !(k in context)) { context[k] = ctx }
    }
    BEGIN { print "chr", "pos", "strand", "context", "site_id", "up_in_WT_27",
                  "trua1b_dependent_27", "combination" }
    ARGIND == 1 {
      if (FNR == 1 || $9 != "up_in_WT_27") { next }
      gsub(/\r/, "")
      k = $1 SUBSEP $2 SUBSEP $3
      register(k, $1, $2, $3, $4)
      up[k] = 1
      next
    }
    ARGIND == 2 {
      if (FNR == 1) { next }
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
  local rna_type="$1"   # mRNA / rRNA
  local up_tsv="${RESPONSE_ROOT}/${rna_type}/WT_27_vs_WT_21.response.tsv"
  local dep_tsv="${BB_ROOT}/${rna_type}/27/trua1b/trua1b-27.dependent.tsv"

  require_file "$up_tsv"
  require_file "$dep_tsv"

  write_overlap_table "$up_tsv" "$dep_tsv" \
    "${SCRIPT_DIR}/${rna_type}_27up_vs_trua1b_dependent_overlap.tsv" "$rna_type"
}

main() {
  write_pair mRNA
  write_pair rRNA
}

main "$@"
