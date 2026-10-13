#!/usr/bin/env bash
# data_prep.sh — 生成本基因型的 tRNA dependent 位点表（组织拆分 + 跨组织并集）
# 职责：只产数据，不画图。逻辑全在共享脚本 ../build_genotype_dependent_tables.sh，
#       本文件只绑定基因型与输出目录。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOTYPE="trua1a"

bash "$SCRIPT_DIR/../build_genotype_dependent_tables.sh" "$GENOTYPE" "$SCRIPT_DIR"

test -s "$SCRIPT_DIR/${GENOTYPE}_all_tissues.dependent.tsv"
test -s "$SCRIPT_DIR/${GENOTYPE}_union.dependent.tsv"
