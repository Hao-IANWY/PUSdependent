#!/usr/bin/env bash
# data_prep.sh — 汇总 polyARNA/STAT/dependent_betabinomial/mRNA 各基因型×组织的 dependent 位点数
# 职责：生成 TSV 数据文件，供 plot.R 读取。不承担画图。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
INPUT_DIR="$ROOT_DIR/polyARNA/STAT/dependent_betabinomial/mRNA"
OUT_TSV="$SCRIPT_DIR/dependent_site_count.tsv"

tissues=("Root" "leaf" "Flower")
genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")

printf 'genotype\ttissue\tcount\n' > "$OUT_TSV"

for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    dep_file="$INPUT_DIR/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
    if [[ -s "$dep_file" ]]; then
      # 每行一个 dependent 位点；首行是表头，减去 1 即为位点数
      count=$(( $(wc -l < "$dep_file") - 1 ))
    else
      count=0
    fi
    printf '%s\t%s\t%s\n' "$genotype" "$tissue" "$count" >> "$OUT_TSV"
  done
done

echo "[data_prep] Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) data rows)"
