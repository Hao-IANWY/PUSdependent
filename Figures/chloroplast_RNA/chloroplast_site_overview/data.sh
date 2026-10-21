#!/usr/bin/env bash
# 生成 WT leaf 叶绿体 rRNA/tRNA/mRNA Ψ 位点总览图。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
