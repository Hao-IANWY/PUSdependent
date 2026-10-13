#!/usr/bin/env bash
# data_prep.sh — 汇总 tRNA/STAT/dependent/ncRNA（纯效应量口径）各基因型×组织的 dependent 位点数
# 职责：生成 TSV 数据文件，供 plot.R 读取。不承担画图。
#
# 口径：只统计 tRNA 来源的位点。tRNA 流程的 ncRNA 参考是
# tRNA-unique_sn_snoRNA.fa（共 589 条：220 条 tRNA「tRNA-aa-codon-n」
# + 287 条 snoRNA「AT…-sno」+ 82 条 snRNA「AT…-sn」），dependent 结果里
# 混有 sn/snoRNA 位点，按 chr 前缀筛掉。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
INPUT_DIR="$ROOT_DIR/tRNA/STAT/dependent/ncRNA"
OUT_TSV="$SCRIPT_DIR/dependent_site_count.tsv"

readonly TRNA_CHR_PATTERN='^tRNA-'

tissues=("Root" "leaf" "Flower")
genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")

printf 'genotype\ttissue\tcount\n' > "$OUT_TSV"

for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    dep_file="$INPUT_DIR/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
    if [[ -s "$dep_file" ]]; then
      # 每行一个 dependent 位点，首行是表头；只数 chr 以 tRNA- 开头的行
      count=$(awk -F'\t' -v pat="$TRNA_CHR_PATTERN" 'NR > 1 && $1 ~ pat {n++} END {print n+0}' "$dep_file")
    else
      count=0
    fi
    printf '%s\t%s\t%s\n' "$genotype" "$tissue" "$count" >> "$OUT_TSV"
  done
done

echo "[data_prep] Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) data rows)"
