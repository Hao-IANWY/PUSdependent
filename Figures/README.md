# Figures Directory Overview

This directory contains scripts and data files for generating publication-quality figures from the analysis of RNA modification landscapes.

## Structure

- `QC/`: Quality control metrics and plots for mRNA and tRNA data.
- `chloroplast_RNA/`: Analysis and visualization of chloroplast RNA data.
- `dependent_betabinomial/`: Visualizations based on Beta-binomial statistical modeling of dependent sites.
- `shared/`: Shared plotting utilities and styles.
- `wildtype_landscape/`: Landscapes of RNA modifications in wild-type samples across different tissues.

## Usage

Each sub-directory typically includes:

- `data.sh` / `data_prep.sh`: Shell scripts to prepare the required input data for plotting.
- `plot.R`: R scripts to generate the figures.
- `*.pdf`, `*.png`, `*.svg`: Generated output figures.
- `*.tsv`: Data files used for plotting.

Refer to the specific `README.md` files within each sub-directory for detailed instructions on how to regenerate individual plots.
