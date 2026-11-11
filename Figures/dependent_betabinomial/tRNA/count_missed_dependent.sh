#!/usr/bin/env bash
# count_missed_dependent.sh — Count tRNA writer-dependent Ψ sites whose decrease is
# observable yet rejected by the significance-q gate.
#
# Background: the beta-binomial model with a global theta retains power in tissues with
# unbalanced sequencing depth (Root). This script quantifies the method's own missed-call
# structure across all genotypes × tissues.
#
# Verdicts (the four gates are replayed per tissue on each genotype's cross-tissue union
# of dependent sites):
#   called         is_dependent TRUE in both alleles (passes the final call)
#   missed_q       status=ok in both alleles, delta>=0.05 and beta_coef<0 in both
#                  (decrease observable) yet not called finally — rejected by the q>=0.05
#                  gate (signal without power at unbalanced depth)
#   missed_depth   status != ok in either allele (filtered by the depth gate / outside
#                  that tissue's candidate universe)
#   missed_effect  the remainder (delta<0.05 or beta_coef>=0 in at least one allele — no
#                  decrease to begin with; mostly tissue-specific true sites, not missed
#                  calls)
# Also reports the median candidate-universe depth of WT and of each mutant allele per
# (genotype, tissue) (the depth-imbalance factor) as a mechanistic control.
#
# Usage: bash count_missed_dependent.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
STAT_DIR="$REPO_ROOT/tRNA/STAT/dependent_betabinomial/ncRNA"
ADJUST_DIR="$REPO_ROOT/tRNA/STAT/adjustGap/ncRNA"

OUT_DETAIL="$SCRIPT_DIR/missed_dependent_by_tissue.tsv"
OUT_SUMMARY="$SCRIPT_DIR/missed_dependent_summary.tsv"

readonly -a GENOTYPES=(trua1a trua1b trua2 trua3 trua4 trua5 trua6)
readonly -a TISSUES=(Root leaf Flower)

# Median of a single column of values
median() { sort -g | awk '{ a[NR] = $1 } END { if (NR == 0) print "NA"; else print (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2 }'; }

printf 'genotype\ttissue\tchr\tpos\tfix_pos\tstrand\tpeak_tissue\tTissues\tn_tissues\tverdict\tl1_status\tl1_delta\tl1_beta\tl1_q\tl2_status\tl2_delta\tl2_beta\tl2_q\n' > "$OUT_DETAIL"
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

    # ---- Gate replay: non-called sites go to OUT_DETAIL, counts go to tmp ----
    tmp_cnt="$(mktemp)"; trap 'rm -f "$tmp_cnt"' EXIT
    awk -F'\t' -v OFS='\t' -v g="$g" -v t="$t" -v detail="$OUT_DETAIL" '
      FNR == 1 { next }                          # Skip the header of each of the three inputs
      FILENAME == ARGV[1] {                      # Union table (tRNA-only)
          key = $1 "|" $2 "|" $4
          u_fix[key] = $3; u_peak[key] = $12; u_tis[key] = $10; u_ntis[key] = $11
          u_chr[key] = $1; u_pos[key] = $2; u_strand[key] = $4
          order[++n_u] = key
          next
      }
      FILENAME == ARGV[2] { k = $1 "|" $2 "|" $3; s1[k]=$14; d1[k]=$8; m1[k]=$11; q1[k]=$10; p1[k]=$15; next }
      FILENAME == ARGV[3] { k = $1 "|" $2 "|" $3; s2[k]=$14; d2[k]=$8; m2[k]=$11; q2[k]=$10; p2[k]=$15; next }
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
          else if (dd1 >= 0.05 && dd2 >= 0.05 && mm1 < 0 && mm2 < 0)
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

    # ---- Median candidate-universe depth: WT 3 replicates (Treat depth columns) + 3 replicates per allele ----
    # Candidate-universe keys (chr/pos/strand) come from the allele-level dependent-sites.tsv, for grepping adjustGap
    l1_keys="$(mktemp)"; l2_keys="$(mktemp)"
    tail -n +2 "$line1_tsv" | cut -f1-3 | tr -d '\r' | sort -u > "$l1_keys"
    tail -n +2 "$line2_tsv" | cut -f1-3 | tr -d '\r' | sort -u > "$l2_keys"
    wt_vals=(); l1_vals=(); l2_vals=()
    for r in 1 2 3; do
      wt_vals+=("$(awk -F'\t' -v rep="rep${r}" '
        FNR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/\r$/, "", h); if (h ~ "^Treat-.*\\." rep "_depth$") c = i } next }
        { print $(c) + 0 }' "$wt_tsv" | median)")
      l1_vals+=("$(LC_ALL=C grep -F -f "$l1_keys" "$ADJUST_DIR/Treat-${g}_1-${t}.rep${r}.adjusted.tsv" | cut -f6 | median)")
      l2_vals+=("$(LC_ALL=C grep -F -f "$l2_keys" "$ADJUST_DIR/Treat-${g}_2-${t}.rep${r}.adjusted.tsv" | cut -f6 | median)")
    done
    rm -f "$l1_keys" "$l2_keys"
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
echo "[count] summary: genotypes=${#GENOTYPES[@]} total union sites=$total_union"
echo "[count] missed_q (genotype×tissue×site combinations) = $total_missed_q; distinct sites = $n_sites; of which Root = $n_root"
