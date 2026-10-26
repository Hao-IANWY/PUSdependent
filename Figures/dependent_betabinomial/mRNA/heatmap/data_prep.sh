#!/usr/bin/env bash
# data_prep.sh — 筛选 WT 修饰比例 >= 20% 的 mRNA dependent 位点并汇总跨基因型比例
#
# 位点集合：每个组织内，7 个突变体的最终 dependent.tsv 并集；WT_Ratio >= 0.20。
# 比例口径：WT 取 dependent 结果自带的 WT_Ratio；突变体取两条株系 dependent-sites.tsv
#           中 Mut_Ratio 的算术均值，与 dependent 鉴定及 loss_effect 图组一致。
# 输出：
#   dependent_sites_wt20.tsv        入选位点及其实际 dependent 基因型（每组织一行/位点）
#   dependent_site_ratio_heatmap.tsv  位点 × (WT + 7 突变体) 的长表，供 plot.R 读取
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
INPUT_DIR="$ROOT_DIR/polyARNA/STAT/dependent_betabinomial/mRNA"
SELECTED_TSV="$SCRIPT_DIR/dependent_sites_wt20.tsv"
OUT_TSV="$SCRIPT_DIR/dependent_site_ratio_heatmap.tsv"
WT_THRESHOLD=0.20

genotypes=(trua1a trua1b trua2 trua3 trua4 trua5 trua6)
tissues=(Root leaf Flower)

tmpdir="$(mktemp -d)"
trap 'rm -rf "$tmpdir"' EXIT
calls_tsv="$tmpdir/dependent_calls.tsv"
: > "$calls_tsv"

# 最终 dependent 表提供位点并集、统一的 WT 比例和实际 dependent 基因型。
for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    dep_file="$INPUT_DIR/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
    [[ -s "$dep_file" ]] || {
      echo "[data_prep] Missing or empty input: $dep_file" >&2
      exit 1
    }
    awk -F '\t' -v OFS='\t' -v tissue="$tissue" -v genotype="$genotype" '
      NR == 1 {
        for (i = 1; i <= NF; i++) col[$i] = i
        required = "chr pos strand context WT_Ratio"
        n = split(required, req, " ")
        for (i = 1; i <= n; i++) {
          if (!(req[i] in col)) {
            printf("Required column missing in %s: %s\n", FILENAME, req[i]) > "/dev/stderr"
            exit 1
          }
        }
        next
      }
      {
        print tissue, genotype, $(col["chr"]), $(col["pos"]),
              $(col["strand"]), $(col["context"]), $(col["WT_Ratio"])
      }
    ' "$dep_file" >> "$calls_tsv"
  done
done

# 合并同一组织内不同基因型的 dependent 记录；核对 WT_Ratio 在各表间一致。
printf 'tissue\tsite_id\tchr\tpos\tstrand\tcontext\tWT_Ratio\tdependent_genotypes\n' > "$SELECTED_TSV"
awk -F '\t' -v OFS='\t' -v threshold="$WT_THRESHOLD" '
  BEGIN { tissue_rank["Root"] = 1; tissue_rank["leaf"] = 2; tissue_rank["Flower"] = 3 }
  {
    tissue = $1; genotype = $2; chr = $3; pos = $4; strand = $5
    key = tissue SUBSEP chr SUBSEP pos SUBSEP strand
    wt_now = $7 + 0
    if (!(key in wt)) {
      wt[key] = wt_now
      context[key] = $6
      dependent_genotypes[key] = genotype
    } else {
      if ((wt[key] - wt_now > 1e-8) || (wt_now - wt[key] > 1e-8)) {
        printf("WT_Ratio mismatch for %s:%s:%s:%s\n", tissue, chr, pos, strand) > "/dev/stderr"
        exit 1
      }
      dependent_genotypes[key] = dependent_genotypes[key] "," genotype
    }
    site_id[key] = chr ":" pos ":" strand
    chr_map[key] = chr; pos_map[key] = pos; strand_map[key] = strand
  }
  END {
    for (key in wt) {
      if (wt[key] >= threshold) {
        split(key, fields, SUBSEP)
        tissue = fields[1]
        printf "%d\t%s\t%s\t%s\t%s\t%s\t%s\t%.15g\t%s\n",
               tissue_rank[tissue], tissue, site_id[key], chr_map[key], pos_map[key],
               strand_map[key], context[key], wt[key], dependent_genotypes[key]
      }
    }
  }
' "$calls_tsv" |
  sort -t $'\t' -k1,1n -k8,8gr -k4,4 -k5,5n |
  cut -f2- >> "$SELECTED_TSV"

# 每行对应一个位点在一个样本组中的修饰比例；WT 只占一列，不参与 dependent 标记。
printf 'tissue\tsite_id\tchr\tpos\tstrand\tcontext\tWT_Ratio\tgenotype\tratio\tis_dependent\n' > "$OUT_TSV"
awk -F '\t' -v OFS='\t' '
  NR > 1 { print $1, $2, $3, $4, $5, $6, $7, "WT", $7, "FALSE" }
' "$SELECTED_TSV" >> "$OUT_TSV"

for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    line1_file="$INPUT_DIR/$tissue/$genotype/${genotype}_1-${tissue}.dependent-sites.tsv"
    line2_file="$INPUT_DIR/$tissue/$genotype/${genotype}_2-${tissue}.dependent-sites.tsv"
    for f in "$line1_file" "$line2_file"; do
      [[ -s "$f" ]] || {
        echo "[data_prep] Missing or empty input: $f" >&2
        exit 1
      }
    done

    awk -F '\t' -v OFS='\t' -v tissue="$tissue" -v genotype="$genotype" '
      NR == FNR {
        if (FNR == 1 || $1 != tissue) next
        key = $3 SUBSEP $4 SUBSEP $5
        if (!(key in selected)) {
          selected[key] = 1
          order[++n] = key
        }
        site_id[key] = $2; chr[key] = $3; pos[key] = $4; strand[key] = $5
        context[key] = $6; wt[key] = $7 + 0; dep_genotypes[key] = $8
        next
      }
      FNR == 1 {
        file_index++
        delete col
        for (i = 1; i <= NF; i++) col[$i] = i
        required = "chr pos strand WT_Ratio Mut_Ratio"
        n_req = split(required, req, " ")
        for (i = 1; i <= n_req; i++) {
          if (!(req[i] in col)) {
            printf("Required column missing in %s: %s\n", FILENAME, req[i]) > "/dev/stderr"
            exit 1
          }
        }
        next
      }
      {
        chr_now = $(col["chr"]); pos_now = $(col["pos"]); strand_now = $(col["strand"])
        key = chr_now SUBSEP pos_now SUBSEP strand_now
        if (!(key in selected)) next
        seen_key = key SUBSEP file_index
        if (seen_key in seen) {
          printf("Duplicate site in %s: %s:%s:%s\n", FILENAME, chr_now, pos_now, strand_now) > "/dev/stderr"
          exit 1
        }
        seen[seen_key] = 1
        wt_now = $(col["WT_Ratio"]) + 0
        if ((wt[key] - wt_now > 1e-6) || (wt_now - wt[key] > 1e-6)) {
          printf("WT_Ratio mismatch for %s %s:%s:%s in %s\n",
                 genotype, chr_now, pos_now, strand_now, FILENAME) > "/dev/stderr"
          exit 1
        }
        mut_sum[key] += $(col["Mut_Ratio"]) + 0
        mut_n[key]++
      }
      END {
        for (i = 1; i <= n; i++) {
          key = order[i]
          is_dependent = (index("," dep_genotypes[key] ",", "," genotype ",") > 0 ? "TRUE" : "FALSE")
          if (mut_n[key] == 2) {
            ratio = mut_sum[key] / 2
            if (ratio < 0 || ratio > 1) {
              printf("Ratio outside [0,1] for %s %s\n", genotype, site_id[key]) > "/dev/stderr"
              exit 1
            }
            ratio_out = sprintf("%.15g", ratio)
          } else {
            # 株系记录不完整时保留 NA，绝不把缺失数据解释为 0% 修饰。
            ratio_out = "NA"
            if (is_dependent == "TRUE") {
              printf("Dependent site lacks two mutant-line values: %s %s %s\n",
                     tissue, genotype, site_id[key]) > "/dev/stderr"
              exit 1
            }
            printf("[data_prep] Warning: %s %s %s has %d/2 mutant-line values; ratio=NA\n",
                   tissue, genotype, site_id[key], mut_n[key]) > "/dev/stderr"
          }
          print tissue, site_id[key], chr[key], pos[key], strand[key], context[key],
                sprintf("%.15g", wt[key]), genotype, ratio_out, is_dependent
        }
      }
    ' "$SELECTED_TSV" "$line1_file" "$line2_file" >> "$OUT_TSV"
  done
done

selected_count=$(( $(wc -l < "$SELECTED_TSV") - 1 ))
actual_rows=$(( $(wc -l < "$OUT_TSV") - 1 ))
expected_rows=$(( selected_count * 8 ))
if [[ "$actual_rows" -ne "$expected_rows" ]]; then
  echo "[data_prep] Row-count check failed: expected $expected_rows, got $actual_rows" >&2
  exit 1
fi

for tissue in "${tissues[@]}"; do
  count=$(awk -F '\t' -v tissue="$tissue" 'NR > 1 && $1 == tissue { n++ } END { print n + 0 }' "$SELECTED_TSV")
  echo "[data_prep] $tissue: $count selected sites (WT_Ratio >= $WT_THRESHOLD)"
done
echo "[data_prep] Generated: $SELECTED_TSV ($selected_count tissue-site rows)"
echo "[data_prep] Generated: $OUT_TSV ($actual_rows ratio rows)"
