#!/usr/bin/env Rscript
# ============================================================================
# plot.R — Scatter plots of WT temperature-responsive Ψ sites at 21 / 27 °C
#          (mRNA / rRNA, one each)
#
# Input: <rna>_temperature_responsive_scatter.tsv (built by data_prep.sh)
#   One row per candidate site (union of Fisher-positive sites at both temperatures),
#   carrying the modification fraction at each temperature and the regulation class
#   (non-responsive / down-regulated / up-regulated).
# Output: <rna>_temperature_responsive_scatter.pdf / .png / .svg
#   X axis = WT modification fraction at 21 °C, Y axis = at 27 °C, dashed line = y = x;
#   grey points = non-responsive sites, red = up-regulated at 27 °C, blue =
#   down-regulated at 27 °C.
#
# Styling comes from Figures/shared. The tRNA version lives in the sibling directory
# Temperature_tRNA/.
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

rna_types <- c("mRNA", "rRNA")

ratio_breaks <- seq(0, 1, by = 0.2)
# Both panels share one legend order and palette; drop = FALSE keeps the legend at three
# complete rows even when a class is empty, so the two figures can be compared
# side by side.
regulation_levels <- c("non-responsive", "down-regulated", "up-regulated")
regulation_colors <- c(
  "non-responsive" = "#B8B8B8CC",
  "down-regulated" = "#3C5488FF",
  "up-regulated"   = "#E64B35FF"
)

init_theme(base_size = SIZE_TOKENS$base_size, sizes = SIZE_TOKENS)

build_plot <- function(df) {
  # ggplot2 4.x draws no legend color key for empty factor levels under drop = FALSE
  # (the row would be text only). Seed the missing levels with one NA-coordinate row to
  # generate the legend key; no point is actually drawn (na.rm = TRUE).
  missing_levels <- setdiff(regulation_levels, as.character(unique(df$regulation)))
  legend_seed <- data.frame(
    WT_21_Ratio = rep(NA_real_, length(missing_levels)),
    WT_27_Ratio = rep(NA_real_, length(missing_levels)),
    regulation = factor(missing_levels, levels = regulation_levels)
  )
  background_df <- bind_rows(
    filter(df, .data$regulation == "non-responsive"),
    legend_seed
  )

  # draw the grey points first and the responsive sites last so responsive sites sit on top
  ggplot(df, aes(x = .data$WT_21_Ratio, y = .data$WT_27_Ratio)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = SIZE_TOKENS$base_linewidth,
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_point(
      data = background_df,
      aes(fill = .data$regulation),
      shape = 21, size = 3.0, color = "black", stroke = 0.5, alpha = 0.88,
      na.rm = TRUE
    ) +
    geom_point(
      data = filter(df, .data$regulation != "non-responsive"),
      aes(fill = .data$regulation),
      shape = 21, size = 3.0, color = "black", stroke = 0.5, alpha = 0.88
    ) +
    annotate(
      "text",
      x = 0.03,
      y = 0.97,
      label = "27°C vs. 21°C",
      hjust = 0,
      vjust = 1,
      size = SIZE_TOKENS$annotation_text_size / ggplot2::.pt,
      family = "Arial",
      fontface = "bold",
      color = "black"
    ) +
    scale_fill_manual(
      values = regulation_colors,
      breaks = regulation_levels,
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
      x = expression(bold("WT 21") * degree * bold("C Modification Ratio")),
      y = expression(bold("WT 27") * degree * bold("C Modification Ratio"))
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
  input_file <- file.path(script_dir, paste0(rna_type, "_temperature_responsive_scatter.tsv"))
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
      WT_21_Ratio = as.numeric(WT_21_Ratio),
      WT_27_Ratio = as.numeric(WT_27_Ratio),
      Delta_Ratio_27_minus_21 = as.numeric(Delta_Ratio_27_minus_21),
      regulation = factor(regulation, levels = regulation_levels)
    )

  if (nrow(plot_df) == 0) {
    stop("No rows found in plotting table: ", input_file)
  }
  if (any(is.na(plot_df$regulation))) {
    stop("Unexpected regulation value in: ", input_file)
  }

  n_responsive <- sum(plot_df$regulation != "non-responsive")
  if (n_responsive == 0) {
    warning(rna_type, ": no responsive site found; the figure will be all grey.")
  }

  save_plot(
    build_plot(plot_df),
    file.path(script_dir, paste0(rna_type, "_temperature_responsive_scatter")),
    width = 6.8,
    height = 6.8
  )
  message(sprintf("[plot] %s: %d candidate sites (%d responsive)", rna_type, nrow(plot_df), n_responsive))
}

message("[plot] Done: mRNA / rRNA temperature-responsive scatter plots (pdf/png/svg)")
