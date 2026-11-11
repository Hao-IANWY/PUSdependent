#!/usr/bin/env Rscript
# ============================================================================
# plot.R — Dot plots of the loss of modification at trua1b-dependent Ψ sites
#          (tRNA × 21 °C / 27 °C, 2 panels)
#
# Input: <rna>_<temp>_dependent_site_loss_effect.tsv (built by data_prep.sh)
#   One row per WT-positive site at that temperature: WT_Ratio / Mut_Ratio /
#   site_status
# Output: <rna>_<temp>_dependent_site_loss_effect.pdf / .png / .svg
#   X axis = mean mutant modification fraction, Y axis = WT modification fraction,
#   dashed line = y = x; grey points = all WT-positive sites, red points = trua1b
#   dependent sites at that temperature.
#   A site below the diagonal (WT > Mut) indicates a TRUA1B-dependent loss of
#   modification.
#
# tRNA version (2 panels); the mRNA / rRNA version lives in the sibling directory
# Temperature_mRNA/. Styling comes from Figures/shared.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

rna_types <- "tRNA"
temps <- c("21", "27")

ratio_breaks <- seq(0, 1, by = 0.2)
site_status_levels <- c("wt_only", "dependent")
site_status_labels <- c("WT sites", "Dependent sites")
status_colors <- setNames(PAL_STATUS, site_status_labels)

init_theme(base_size = SIZE_TOKENS$base_size, sizes = SIZE_TOKENS)

build_plot <- function(df, temp_label) {
  ggplot(df, aes(x = .data$Mut_Ratio, y = .data$WT_Ratio)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = SIZE_TOKENS$base_linewidth,
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_point(
      data = filter(df, .data$site_status == "WT sites"),
      aes(fill = .data$site_status),
      shape = 21, size = 2.8, color = "black", stroke = 0.45
    ) +
    geom_point(
      data = filter(df, .data$site_status == "Dependent sites"),
      aes(fill = .data$site_status),
      shape = 21, size = 3.0, color = "black", stroke = 0.55
    ) +
    annotate(
      "text",
      x = 0.03,
      y = 0.97,
      label = temp_label,
      hjust = 0,
      vjust = 1,
      size = SIZE_TOKENS$annotation_text_size / ggplot2::.pt,
      family = "Arial",
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(
      values = status_colors,
      breaks = site_status_labels,
      drop = FALSE,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = ratio_breaks,
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = ratio_breaks,
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = 0.02)
    ) +
    coord_equal() +
    labs(
      x = expression(bolditalic("trua1b") * bold(" Mutant Modification Ratio")),
      y = expression(bold("WT Modification Ratio"))
    ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", colour = NA),
      axis.title = element_text(face = "bold", colour = "black"),
      plot.margin = margin(10, 20, 10, 12)
    )
}

for (rna_type in rna_types) {
  for (temp in temps) {
    stem <- paste0(rna_type, "_", temp, "_dependent_site_loss_effect")
    input_file <- file.path(script_dir, paste0(stem, ".tsv"))
    if (!file.exists(input_file)) {
      stop("Input table not found: ", input_file, "\nRun data_prep.sh first.")
    }

    plot_df <- read.table(
      input_file,
      header = TRUE,
      sep = "\t",
      quote = "",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    ) |>
      mutate(
        WT_Ratio = as.numeric(WT_Ratio),
        Mut_Ratio = as.numeric(Mut_Ratio),
        Delta_Ratio = as.numeric(Delta_Ratio),
        site_status = factor(site_status,
                             levels = site_status_levels,
                             labels = site_status_labels)
      )

    if (nrow(plot_df) == 0) {
      stop("No rows found in plotting table: ", input_file)
    }

    n_dependent <- sum(plot_df$site_status == "Dependent sites")
    if (n_dependent == 0) {
      # tRNA can lack dependent sites at both temperatures, leaving an all-grey figure;
      # render it anyway but log a note so this is not later mistaken for missing data
      message("[plot] ", rna_type, "/", temp, ": no dependent site; the figure is all grey.")
    }

    save_plot(
      build_plot(plot_df, paste0(temp, "°C")),
      file.path(script_dir, stem),
      width = 6.8,
      height = 6.8
    )
    message(sprintf("[plot] %s/%s: %d WT sites (%d dependent)", rna_type, temp, nrow(plot_df), n_dependent))
  }
}

message("[plot] Done: tRNA trua1b dependent-site loss-effect plots (pdf/png/svg)")
