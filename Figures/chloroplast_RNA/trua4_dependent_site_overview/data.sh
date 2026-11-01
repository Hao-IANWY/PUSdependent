#!/usr/bin/env bash
# 生成 Leaf 叶绿体 rRNA / tRNA / mRNA 位点 WT vs trua4 散点图与汇总表。
# 数据源为 STAT/dependent_betabinomial（beta-binomial 口径），dependent 位点
# 采用两株系 AND 判定，与 Figures/dependent_betabinomial 一致。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
