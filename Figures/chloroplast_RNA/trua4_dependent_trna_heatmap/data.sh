#!/usr/bin/env bash
# 生成 trua4 dependent 叶绿体 tRNA 位点的结构热图。
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
