#!/usr/bin/env bash
# data_prep.sh — 汇总 7 个 trua 突变体 × 3 个组织的 mRNA dependent 位点修饰下降数据
# 职责：生成 dependent_site_loss_effect_dotplot.tsv，供 plot.R 读取。不承担画图。
#
# 每行一个 WT 阳性 mRNA 位点（背景全量），标注其是否属于对应基因型×组织的 dependent 集合：
#   WT_Ratio  来自 dependent-sites.tsv 的 WT_Ratio（WT 阳性候选的 Treat 平均修饰比）
#   Mut_Ratio 来自 dependent-sites.tsv 的 Mut_Ratio（adjustGap 3 个 Treat 重复的均值，未筛选），两株系再平均
#
# 注意：Mut_Ratio 必须取未筛选的 adjustGap 口径，而非 fisher 表。fisher 会把突变体中
# 低信号位点过滤掉、置为 0，但这些位点其实有微弱残留信号（如 AT1G04210.1 在 adjustGap
# 中 ratio=0.000/0.009174/0.012658 → 均值 0.007277，而 fisher 里它被过滤、Mut_Ratio=0）。
# dependent-sites.tsv 的 Mut_Ratio 正是 call_dependent.py 从 adjustGap 3 个 Treat 重复算得的均值。
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
DEPENDENT_DIR="$ROOT_DIR/polyARNA/STAT/dependent/mRNA"
OUT_TSV="$SCRIPT_DIR/dependent_site_loss_effect_dotplot.tsv"

genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")
tissues=("Root" "leaf" "Flower")

printf 'panel_id\tgenotype\ttissue\tchr\tpos\tstrand\tcontext\tsite_id\tWT_Ratio\tMut_Ratio\tDelta_Ratio\tsite_status\n' > "$OUT_TSV"

for genotype in "${genotypes[@]}"; do
  for tissue in "${tissues[@]}"; do
    line1_file="$DEPENDENT_DIR/$tissue/$genotype/${genotype}_1-${tissue}.dependent-sites.tsv"
    line2_file="$DEPENDENT_DIR/$tissue/$genotype/${genotype}_2-${tissue}.dependent-sites.tsv"

    for f in "$line1_file" "$line2_file"; do
      [[ -s "$f" ]] || { echo "[data_prep] Missing input file: $f" >&2; exit 1; }
    done

    awk -F '\t' -v OFS='\t' \
        -v genotype="$genotype" \
        -v tissue="$tissue" \
        -v panel_id="${genotype}_${tissue}" '
        function validate_required(required_names,    i, name) {
            split(required_names, fields, ",")
            for (i in fields) {
                name = fields[i]
                if (!(name in col)) {
                    printf("Required column missing in %s: %s\n", FILENAME, name) > "/dev/stderr"
                    exit 1
                }
            }
        }
        FNR == 1 {
            file_index += 1
            delete col
            for (i = 1; i <= NF; i++) {
                sub(/\r$/, "", $i)
                col[$i] = i
            }
            validate_required("chr,pos,strand,context,WT_Ratio,Mut_Ratio,is_dependent")
            next
        }
        {
            chr = $(col["chr"])
            pos = $(col["pos"])
            strand = $(col["strand"])
            site_id = chr "|" pos "|" strand
            wt_ratio = $(col["WT_Ratio"])
            mut_ratio = $(col["Mut_Ratio"])
            dep_raw = $(col["is_dependent"])
            sub(/\r$/, "", dep_raw)
            is_dep = (dep_raw == "True" ? 1 : 0)

            if (file_index == 1) {
                n += 1
                order[n] = site_id
                chr_map[site_id] = chr
                pos_map[site_id] = pos
                strand_map[site_id] = strand
                context_map[site_id] = $(col["context"])
                wt_map[site_id] = wt_ratio
                mut_sum[site_id] = mut_ratio
                mut_count[site_id] = 1
                dep_map[site_id] = is_dep
                next
            }

            # file_index == 2：与 line1 合并；累加 Mut_Ratio、两线 is_dependent 同时成立才算 dependent
            if (site_id in wt_map) {
                mut_sum[site_id] += mut_ratio
                mut_count[site_id] += 1
                dep_map[site_id] = (dep_map[site_id] && is_dep)
            }
        }
        END {
            for (i = 1; i <= n; i++) {
                site_id = order[i]
                mut_ratio = (mut_count[site_id] > 0 ? mut_sum[site_id] / mut_count[site_id] : 0)
                delta = wt_map[site_id] - mut_ratio
                status = (dep_map[site_id] ? "dependent" : "wt_only")
                print panel_id, genotype, tissue, chr_map[site_id], pos_map[site_id], strand_map[site_id], \
                      context_map[site_id], site_id, wt_map[site_id], mut_ratio, delta, status
            }
        }
    ' "$line1_file" "$line2_file" >> "$OUT_TSV"

    echo "[data_prep] $genotype / $tissue done"
  done
done

echo "[data_prep] Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) data rows)"
