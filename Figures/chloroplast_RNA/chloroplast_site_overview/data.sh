#!/usr/bin/env bash
# Render the overview of chloroplast rRNA/tRNA/mRNA Ψ sites in WT leaf.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
