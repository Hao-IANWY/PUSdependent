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

genotype <- "trua2"

# ---- Tissue-split heatmap: rows = tissue × representative tRNA, value = Delta Ratio (WT - Mut). Produced for all 7 genotypes ----
tissue_result <- run_tissue_split_heatmap(
  genotype = genotype,
  base_dir = tRNA_dir,
  out_dir = script_dir
)
print(tissue_result$metadata)

source(file.path(tRNA_dir, "dependent_trna_distribution_common.R"))

# ---- Ψ site distribution box plot: x = tissue, y = modification fraction (%), three groups WT / <genotype>_1 / <genotype>_2 ----
distribution_result <- run_single_panel_dependent_position_distributions(
  genotype = genotype,
  input_path = file.path(script_dir, paste0(genotype, "_all_tissues.dependent.tsv")),
  out_dir = script_dir,
  plot_specs = tibble::tibble(
    compartment = c("Mt"),
    label = c(paste0("mt-tRNA \u03a8", "27/28")),
    fix_pos_values = list(c(27L, 28L)),
    output_suffix = c("mt_27_28")
  )
)
print(distribution_result$metadata)

source(file.path(tRNA_dir, "dependent_trna_motif_common.R"))

# ---- Site motif seqlogo: 5-nt contexts of union-deduplicated sites, selected by compartment × Sprinzl position ----
motif_result <- run_single_panel_dependent_motifs(
  genotype = genotype,
  input_path = file.path(script_dir, paste0(genotype, "_union.dependent.tsv")),
  out_dir = script_dir,
  plot_specs = tibble::tibble(
    compartment = c("Mt"),
    label = c(paste0("mt-tRNA Ψ", "27/28")),
    fix_pos_values = list(c(27L, 28L)),
    output_suffix = c("mt_27_28")
  )
)
print(motif_result$metadata)
