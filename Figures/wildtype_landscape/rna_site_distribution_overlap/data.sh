#!/usr/bin/env bash
# Render the tissue-wise Ψ-site density panels and the pairwise site-density
# overlap correlation figure via plot.R.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
