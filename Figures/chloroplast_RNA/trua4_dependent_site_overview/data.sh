#!/usr/bin/env bash
# 生成叶绿体 WT vs trua4 位点散点图（rRNA / tRNA / mRNA）。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
