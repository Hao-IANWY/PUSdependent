#!/usr/bin/env Rscript
# plot.R — Sprinzl position 28 figure set (Flower panel; one heatmap + one
# modification-fraction distribution)
#
# Input (built by data_prep.sh): pos28_flower_site_table.tsv
# The implementation lives in ../site_specificity_common.R; this file only declares
# parameters.
# Panels display "Flowers" following the project-wide display convention; the
# underlying directories and data keep Flower.

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
common_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
shared_dir <- normalizePath(file.path(common_dir, "..", "..", "..", "shared"), mustWork = TRUE)

source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))
source(file.path(common_dir, "site_specificity_common.R"))

genotypes <- c("trua1a", "trua5")
sprinzl_pos <- 28L
position_label <- paste0("Ψ", sprinzl_pos)
input_path <- file.path(script_dir, "pos28_flower_site_table.tsv")
file_prefix <- "pos28_flower"

# ---- Heatmap: rows = nuclear representative tRNAs modified at Ψ28 in WT (descending
# WT fraction), columns = WT + each allele ----
heatmap_result <- run_site_specificity_heatmap(
  input_path = input_path,
  genotypes = genotypes,
  out_dir = script_dir,
  file_prefix = file_prefix
)
print(heatmap_result$metadata)

# ---- Distribution: x = WT / trua1a_1 / trua1a_2 / trua5_1 / trua5_2 ----
distribution_result <- run_site_specificity_distribution(
  input_path = input_path,
  genotypes = genotypes,
  position_label = position_label,
  out_dir = script_dir,
  file_prefix = file_prefix
)
print(distribution_result$summary)
