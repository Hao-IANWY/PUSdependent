#!/usr/bin/env bash
# Render the circular overview of Ψ sites on the chloroplast chromosome
# (plain version plus an mRNA-gene-annotated version).
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
