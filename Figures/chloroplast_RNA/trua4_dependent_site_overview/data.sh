#!/usr/bin/env bash
# Render WT vs trua4 scatter plots and summary tables for chloroplast
# rRNA / tRNA / mRNA sites in leaf.
# Data come from STAT/dependent_betabinomial (beta-binomial calling criteria);
# dependent sites follow the AND rule across both mutant alleles, matching
# Figures/dependent_betabinomial.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
