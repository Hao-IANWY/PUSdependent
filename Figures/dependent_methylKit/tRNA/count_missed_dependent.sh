#!/usr/bin/env bash
# count_missed_dependent.sh — 统计「下降看得到、却被显著性 q 卡掉」的 tRNA dependent 位点
#
# 背景：Root 文库 WT 覆盖度比突变株低几十倍，methylKit 加权 GLM 的过离散校正吃掉
# 全部功效，修饰明显下降的位点也判不上 dependent（详见 trua3/root_depth/ 与
# trua3 的 .full. 分布图）。本脚本对全部基因型 × 组织量化这一现象。
#
# 口径（对每个基因型的跨组织并集 dependent 位点，在每个组织里重演四道门）：
#   called         两株系 is_dependent 均为 TRUE（final 判上）
#   missed_q       两株系 status=ok、delta≥0.05、meth.diff>0 全过（下降看得到），
#                  但 final 没判上 —— 被-q≥0.05 卡掉（深度失衡无功效的信号）
#   missed_depth   任一株系 status != ok（深度门被滤 / 不在该组织候选宇宙）
#   missed_effect  其余（至少一个株系 delta<0.05 或 meth.diff≤0 —— 本就无下降，
#                  多为组织特异真位点，不算漏检）
# 另输出每个 (基因型,组织) 的 WT/突变株候选宇宙中位深度（失衡倍数），作为机制对照。
#
# Usage: bash count_missed_dependent.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
STAT_DIR="$REPO_ROOT/tRNA/STAT/dependent_methylKit/ncRNA"

OUT_DETAIL="$SCRIPT_DIR/missed_dependent_by_tissue.tsv"
OUT_SUMMARY="$SCRIPT_DIR/missed_dependent_summary.tsv"

readonly -a GENOTYPES=(trua1a trua1b trua2 trua3 trua4 trua5 trua6)
readonly -a TISSUES=(Root leaf Flower)

# 单列数值的中位数
median() { sort -g | awk '{ a[NR] = $1 } END { if (NR == 0) print "NA"; else print (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2 }'; }

printf 'genotype\ttissue\tchr\tpos\tfix_pos\tstrand\tpeak_tissue\tTissues\tn_tissues\tverdict\tl1_status\tl1_delta\tl1_methdiff\tl1_q\tl2_status\tl2_delta\tl2_methdiff\tl2_q\n' > "$OUT_DETAIL"
printf 'genotype\ttissue\tn_union\tcalled\tmissed_q\tmissed_depth\tmissed_effect\twt_median_depth\tl1_median_depth\tl2_median_depth\timbalance_x_l1\timbalance_x_l2\n' > "$OUT_SUMMARY"

total_missed_q=0
total_union=0

for g in "${GENOTYPES[@]}"; do
  union_tsv="$SCRIPT_DIR/$g/${g}_union.dependent.tsv"
  [[ -s "$union_tsv" ]] || { echo "Missing union table: $union_tsv" >&2; exit 1; }
  n_union=$(( $(wc -l < "$union_tsv") - 1 ))
  total_union=$(( total_union + n_union ))

  for t in "${TISSUES[@]}"; do
    line1_tsv="$STAT_DIR/$t/$g/${g}_1-$t.dependent-sites.tsv"
    line2_tsv="$STAT_DIR/$t/$g/${g}_2-$t.dependent-sites.tsv"
    wt_tsv="$STAT_DIR/$t/$g/WT-${t}_sites.Treat_avg_ratio.tsv"
    [[ -s "$line1_tsv" && -s "$line2_tsv" && -s "$wt_tsv" ]] || { echo "Missing input for $g/$t" >&2; exit 1; }

    # ---- 门重演：明细（未判上的位点）落 OUT_DETAIL，计数落 tmp ----
    tmp_cnt="$(mktemp)"; trap 'rm -f "$tmp_cnt"' EXIT
    awk -F'\t' -v OFS='\t' -v g="$g" -v t="$t" -v detail="$OUT_DETAIL" '
      FNR == 1 { next }                          # 三个输入的表头一律跳过
      FILENAME == ARGV[1] {                      # 并集表（tRNA-only）
          key = $1 "|" $2 "|" $4
          u_fix[key] = $3; u_peak[key] = $12; u_tis[key] = $10; u_ntis[key] = $11
          u_chr[key] = $1; u_pos[key] = $2; u_strand[key] = $4
          order[++n_u] = key
          next
      }
      FILENAME == ARGV[2] { k = $1 "|" $2 "|" $3; s1[k]=$12; d1[k]=$8; m1[k]=$11; q1[k]=$10; p1[k]=$13; next }
      FILENAME == ARGV[3] { k = $1 "|" $2 "|" $3; s2[k]=$12; d2[k]=$8; m2[k]=$11; q2[k]=$10; p2[k]=$13; next }
      END {
        n_called = n_q = n_depth = n_eff = 0
        for (i = 1; i <= n_u; i++) {
          key = order[i]
          st1 = (key in s1) ? s1[key] : "absent"; st2 = (key in s2) ? s2[key] : "absent"
          dd1 = (key in d1) ? d1[key] + 0 : -1;  dd2 = (key in d2) ? d2[key] + 0 : -1
          mm1 = (key in m1) ? m1[key] + 0 : -1;  mm2 = (key in m2) ? m2[key] + 0 : -1
          dep = (p1[key] == "TRUE" && p2[key] == "TRUE")

          if (dep)                              verdict = "called"
          else if (st1 != "ok" || st2 != "ok")  verdict = "missed_depth"
          else if (dd1 >= 0.05 && dd2 >= 0.05 && mm1 > 0 && mm2 > 0)
                                                verdict = "missed_q"
          else                                  verdict = "missed_effect"

          if (verdict == "called") { n_called++; continue }
          if (verdict == "missed_q") n_q++
          else if (verdict == "missed_depth") n_depth++
          else n_eff++

          print g, t, u_chr[key], u_pos[key], u_fix[key], u_strand[key], \
                u_peak[key], u_tis[key], u_ntis[key], verdict, \
                st1, dd1, mm1, (key in q1) ? q1[key] : "NA", \
                st2, dd2, mm2, (key in q2) ? q2[key] : "NA" >> detail
        }
        print n_called, n_q, n_depth, n_eff
      }' "$union_tsv" "$line1_tsv" "$line2_tsv" > "$tmp_cnt"
    read -r n_called n_q n_depth n_eff < "$tmp_cnt"
    rm -f "$tmp_cnt"

    # ---- 候选宇宙中位深度：WT 3 重复（Treat 深度列）+ 两株系各 3 重复 ----
    wt_vals=(); l1_vals=(); l2_vals=()
    for r in 1 2 3; do
      wt_vals+=("$(awk -F'\t' -v rep="rep${r}" '
        FNR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/\r$/, "", h); if (h ~ "^Treat-.*\\." rep "_depth$") c = i } next }
        { print $(c) + 0 }' "$wt_tsv" | median)")
      l1_vals+=("$(cut -f5 "$STAT_DIR/$t/$g/${g}_1-$t.methylkit.samples/${g}_1-$t.rep${r}.methyl.tsv" | median)")
      l2_vals+=("$(cut -f5 "$STAT_DIR/$t/$g/${g}_2-$t.methylkit.samples/${g}_2-$t.rep${r}.methyl.tsv" | median)")
    done
    mean3() { awk -v a="${1}" -v b="${2}" -v c="${3}" 'BEGIN { print (a + b + c) / 3 }'; }
    wt_med="$(mean3 "${wt_vals[0]}" "${wt_vals[1]}" "${wt_vals[2]}")"
    l1_med="$(mean3 "${l1_vals[0]}" "${l1_vals[1]}" "${l1_vals[2]}")"
    l2_med="$(mean3 "${l2_vals[0]}" "${l2_vals[1]}" "${l2_vals[2]}")"

    awk -v g="$g" -v t="$t" -v n="$n_union" \
        -v c="$n_called" -v q="$n_q" -v d="$n_depth" -v e="$n_eff" \
        -v w="$wt_med" -v a="$l1_med" -v b="$l2_med" 'BEGIN {
      printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%.0f\t%.0f\t%.0f\t%.1f\t%.1f\n",
             g, t, n, c, q, d, e, w, a, b, a / w, b / w
    }' >> "$OUT_SUMMARY"

    total_missed_q=$(( total_missed_q + n_q ))
    echo "[count] $g/$t: union=$n_union called=$n_called missed_q=$n_q missed_depth=$n_depth missed_effect=$n_eff"
  done
done

n_root="$(awk -F'\t' 'NR>1 && $10=="missed_q" && $2=="Root"' "$OUT_DETAIL" | wc -l)"
n_sites="$(awk -F'\t' 'NR>1 && $10=="missed_q" { print $1 "|" $3 "|" $4 "|" $6 }' "$OUT_DETAIL" | sort -u | wc -l)"
echo "[count] 汇总: 基因型数=${#GENOTYPES[@]} 并集位点总数=$total_union"
echo "[count] missed_q (基因型×组织×位点组合) = $total_missed_q；去重位点数 = $n_sites；其中 Root 占 $n_root"
