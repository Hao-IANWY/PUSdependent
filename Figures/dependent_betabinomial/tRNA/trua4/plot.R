#!/usr/bin/env Rscript
# plot.R — tRNA TRUA-dependent Ψ site figure set for this genotype
#
# Inputs (built by data_prep.sh):
#   <genotype>_all_tissues.dependent.tsv  per-tissue dependent sites (with fix_pos and each allele's modification fraction)
#   <genotype>_union.dependent.tsv        cross-tissue union (used only by the union heatmap)
# The implementation lives in ../dependent_trna_heatmap_common.R and
# ../dependent_trna_distribution_common.R; this file only declares the genotype and
# selects which figures to produce.

suppressPackageStartupMessages({
  library(tibble)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
tRNA_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
shared_dir <- normalizePath(file.path(tRNA_dir, "..", "..", "shared"), mustWork = TRUE)

source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))
source(file.path(tRNA_dir, "dependent_trna_heatmap_common.R"))

genotype <- "trua4"

# ---- Tissue-split heatmap: rows = tissue × representative tRNA, value = Delta Ratio (WT - Mut). Produced for all 7 genotypes ----
tissue_result <- run_tissue_split_heatmap(
  genotype = genotype,
  base_dir = tRNA_dir,
  out_dir = script_dir
)
print(tissue_result$metadata)

# ---- Cross-tissue union heatmap: rows = representative tRNA, value = the largest cross-tissue Delta Ratio at that site ----
union_result <- run_union_dependent_heatmap(
  genotype = genotype,
  base_dir = tRNA_dir,
  out_dir = script_dir
)
print(union_result$metadata)
