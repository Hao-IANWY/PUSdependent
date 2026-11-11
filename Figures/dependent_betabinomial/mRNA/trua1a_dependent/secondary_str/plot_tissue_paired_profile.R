#!/usr/bin/env Rscript
# plot_tissue_paired_profile.R — Per-tissue versions of the extended −20 to +60 nt
# paired-probability window for the TRUA1A-dependent Ψ sites at 21 °C, one plot per
# tissue (Root/Leaf/Flower) in the project tissue color, mirroring the single-group
# styling of plot_union_paired_profile.R.
# Input: trua1a_groups_structure_unpaired_profile.tsv (run data_prep.sh first)
# Output: window_-20_to_60/trua1a_<tissue>_-20_to_60_paired.pdf/png/svg
suppressPackageStartupMessages(library(ggplot2))

# Keep Unicode labels intact in Cairo outputs when the shell starts in C.UTF-8.
if (!nzchar(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))) {
  stop("Could not enable an installed UTF-8 locale (en_US.UTF-8).")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

profile <- read.table(file.path(script_dir, "trua1a_groups_structure_unpaired_profile.tsv"),
                      header = TRUE, sep = "\t", quote = "", comment.char = "",
                      stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
if (nrow(profile) == 0) stop("Empty RNAfold input table.")

# Per-tissue sets only; the Union figure lives in plot_union_paired_profile.R.
tissue_sets <- c(Root = "Root", leaf = "Leaf", Flower = "Flower")

sem <- function(x) if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else 0

reg_dir <- file.path(script_dir, "window_-20_to_60")
dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)

for (tissue in names(tissue_sets)) {
  set_name <- unname(tissue_sets[[tissue]])
  sub <- profile[profile$set == set_name & profile$window == 100 &
                 as.numeric(profile$temp_C) == 21, ]
  if (nrow(sub) == 0) stop("No ±100-nt rows found at 21 °C for set: ", set_name)
  n_sites <- nrow(unique(sub[, c("transcript", "pos")]))
  sub$unpaired_prob <- as.numeric(sub$unpaired_prob)
  sub$rel_pos <- as.numeric(sub$rel_pos)
  sub$paired_prob <- 1 - sub$unpaired_prob

  prof_mean <- aggregate(paired_prob ~ rel_pos, data = sub, FUN = mean)
  prof_sem <- aggregate(paired_prob ~ rel_pos, data = sub, FUN = sem)
  names(prof_mean)[2] <- "mean"
  names(prof_sem)[2] <- "sem"
  prof <- merge(prof_mean, prof_sem, by = "rel_pos", sort = FALSE)

  col_tissue <- PAL_TISSUE_NPG[[set_name]]

  init_theme(base_size = 16)

  p <- ggplot(prof[prof$rel_pos >= -20 & prof$rel_pos <= 60, ],
              aes(x = rel_pos, y = mean)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, fill = col_tissue) +
    geom_line(linewidth = 0.75, color = col_tissue) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 1), xlim = c(-20, 60)) +
    scale_x_continuous(breaks = seq(-20, 60, 10)) +
    labs(x = "Position relative to Ψ site (nt)", y = "Mean paired probability",
         title = sprintf("TRUA1A %s-dependent Ψ sites (n=%d)", set_name, n_sites)) +
    theme(plot.title = element_text(size = 12))
  save_plot(p, file.path(reg_dir, paste0("trua1a_", tissue, "_-20_to_60_paired")),
            width = 9.5, height = 5.0)
  message(sprintf("[plot] %s: n=%d -> trua1a_%s_-20_to_60_paired", set_name, n_sites, tissue))
}

message("[plot] Done: per-tissue -20 to +60 nt paired-probability profiles")
