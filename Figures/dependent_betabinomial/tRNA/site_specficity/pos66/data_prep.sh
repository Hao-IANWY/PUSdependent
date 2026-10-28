#!/usr/bin/env bash
# data_prep.sh — Sprinzl 66 位点表（Flower 组）
# 问题：WT Flower 里在 Sprinzl 66 位上有修饰的代表（细胞核来源）tRNA 有多少个？
#       这些 tRNA 在 trua1a 两条株系里的修饰比例如何？
# 分室限定 Nu（老板 2026-10-28 拍板：忽略 Mt/Pt 只保留核）；
# 逻辑全在共享脚本 ../build_site_specificity_table.sh，本文件只绑定参数。
# 加组织：复制本目录改 TISSUE 与文件名；加位号：复制本目录改 POS。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON_DIR="$(dirname "$SCRIPT_DIR")"

TISSUE="Flower"
GENOTYPES="trua1a"
POS=66
COMPARTMENTS="Nu"

bash "$COMMON_DIR/build_site_specificity_table.sh" "$TISSUE" "$GENOTYPES" "$POS" "$SCRIPT_DIR" "$COMPARTMENTS"

test -s "$SCRIPT_DIR/pos${POS}_flower_site_table.tsv"
