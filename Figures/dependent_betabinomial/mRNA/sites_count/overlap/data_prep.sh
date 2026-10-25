#!/usr/bin/env bash
# data_prep.sh — 汇总两种方法 (dependent_betabinomial vs dependent 纯效应量版) 的 dependent 位点集合
# 职责：生成 TSV 数据文件，供 plot.R 读取。不承担画图。
# 位点键：chr:pos:strand（每个 overlap 目录只对应一种 RNA 类型，命名空间一致）
# 粒度：tissue × genotype × method 的去重位点集（供每组织×基因型的 Venn 统计）
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# overlap -> sites_count -> <RNA> -> dependent_betabinomial -> Figures -> 仓库根
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
RNA_TYPE="$(basename "$(cd "$SCRIPT_DIR/../.." && pwd)")"   # mRNA | rRNA

BETABINOMIAL_DIR="$ROOT_DIR/polyARNA/STAT/dependent_betabinomial/$RNA_TYPE"
DEPENDENT_DIR="$ROOT_DIR/polyARNA/STAT/dependent/$RNA_TYPE"
OUT_TSV="$SCRIPT_DIR/dependent_site_overlap.tsv"

tissues=("Root" "leaf" "Flower")
genotypes=("trua1a" "trua1b" "trua2" "trua3" "trua4" "trua5" "trua6")

# 提取某文件的位点键（去表头、去 CR、取 chr:pos:strand）
site_keys() {
  local file="$1"
  tail -n +2 "$file" | tr -d '\r' | awk -F'\t' 'NF>=3 { print $1":"$2":"$3 }'
}

printf 'tissue\tgenotype\tmethod\tsite\n' > "$OUT_TSV"

for tissue in "${tissues[@]}"; do
  for genotype in "${genotypes[@]}"; do
    for method_dir_name in betabinomial dependent; do
      if [[ "$method_dir_name" == "betabinomial" ]]; then
        base="$BETABINOMIAL_DIR"
      else
        base="$DEPENDENT_DIR"
      fi
      dep_file="$base/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
      [[ -s "$dep_file" ]] || continue
      site_keys "$dep_file" | sort -u \
        | awk -v t="$tissue" -v g="$genotype" -v m="$method_dir_name" \
            '{ print t"\t"g"\t"m"\t"$0 }' >> "$OUT_TSV"
    done
  done
done

echo "[data_prep] Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) site rows)"
