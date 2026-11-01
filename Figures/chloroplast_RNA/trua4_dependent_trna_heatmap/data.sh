#!/usr/bin/env bash
# 生成 Leaf 叶绿体 trua4 dependent tRNA 位点的结构热图。
# dependent 位点取 STAT/dependent_betabinomial 的官方 dependent 表（两株系 AND），
# 与 Figures/dependent_betabinomial 口径一致。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
