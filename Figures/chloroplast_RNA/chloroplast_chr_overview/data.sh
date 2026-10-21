#!/usr/bin/env bash
# 生成叶绿体染色体环形 Ψ 位点总览图（普通版 + mRNA 基因标注版）。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
