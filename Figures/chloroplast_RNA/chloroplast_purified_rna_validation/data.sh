#!/usr/bin/env bash
# Render the Ψ-site count bar plot and the poly(A)+ / small RNA overlap Venn
# diagrams for the purified-chloroplast-RNA library.
# Depends on output tables from the sibling chloroplast_site_overview step;
# run that directory's data.sh first.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
