#!/usr/bin/env bash
# Render the structure heatmap of trua4-dependent chloroplast tRNA sites in
# leaf.
# Dependent sites come from the official dependent table of
# STAT/dependent_betabinomial (AND across both mutant alleles), matching the
# calling criteria of Figures/dependent_betabinomial.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
