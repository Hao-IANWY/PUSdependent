#!/usr/bin/env bash
# 生成叶绿体纯化 RNA 位点数量 barplot 与 poly(A)+ / small RNA overlap Venn。
# 依赖同级 chloroplast_site_overview 的产出表，请先运行该目录的 data.sh。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
