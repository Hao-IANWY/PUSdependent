#!/usr/bin/env Rscript

# mRNA metagene site-density figure (metaPlotR-style layout):
# WT positive sites from the three tissues (Root/Leaf/Flower) are projected
# onto normalized 5'UTR-CDS-3'UTR metagene coordinates (written by
# build_metagene_positions.sh; fixed region fractions 15%/70%/15%) and drawn
# as kernel density curves. Dashed lines mark the region boundaries; below
# the panel sits a transcript schematic (thin UTR rectangles + a thick CDS
# rectangle, vertically centered) with region labels; site counts are
# annotated at the top inside the dashed boundary lines; the x axis ends in
# an arrow on the right. The three panels share one y axis (ticks on the
# leftmost panel only); the facet strips keep their text (plural names
# Roots/Leaves/Flowers) but have no frame or fill; no legend.
#
# Inputs:  meta/mrna_metagene_sites.tsv (output of build_metagene_positions.sh)
# Outputs: meta/mrna_metagene_density.{pdf,png,svg}
#          meta/mrna_metagene_density_curve.tsv (density curve data as plotted)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggprism)
  library(tibble)
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]

  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) {
      return(dirname(normalizePath(active_path)))
    }
  }

  normalizePath(getwd())
}

run_metagene_density <- function(out_dir = NULL, output_stem = "mrna_metagene_density", save_outputs = TRUE) {
  script_dir <- get_script_dir()
  # meta sits one level below Figures:
  # meta -> mrna_region_distribution -> wildtype_landscape -> Figures -> repo
  repo_root <- normalizePath(file.path(script_dir, "../../../.."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = 16, border = TRUE, sizes = SIZE_TOKENS)
  out_dir <- if (is.null(out_dir)) script_dir else normalizePath(out_dir, mustWork = FALSE)

  sites_file <- file.path(script_dir, "mrna_metagene_sites.tsv")
  sites <- readr::read_tsv(sites_file, show_col_types = FALSE) |>
    mutate(Tissue = factor(Tissue, levels = c("Root", "Leaf", "Flower")))
  tissue_levels <- levels(sites$Tissue)
  tissue_colors <- setNames(PAL_TISSUE, c("Root", "Leaf", "Flower"))[tissue_levels]

  # Positions of the region-boundary dashed lines; must match U5_FRAC and
  # CDS_FRAC of build_metagene_positions.sh
  region_bounds <- c(0.15, 0.85)

  # Per-tissue kernel density curves (truncated to 0-1, no spill beyond the
  # metagene); one shared y axis whose top tick covers the global peak; the
  # y placement of the schematic and annotations likewise uses fractions of
  # the global peak (a single coordinate system across all three panels)
  curves <- sites |>
    group_by(Tissue) |>
    summarise(fit = list(stats::density(x, from = 0, to = 1, n = 512, cut = 0)),
              n_sites = dplyr::n(), .groups = "drop")

  curve_df <- curves |>
    transmute(Tissue, x = purrr::map(fit, "x"), y = purrr::map(fit, "y")) |>
    unnest(c(x, y))

  ymax <- max(curve_df$y)
  ax <- y_axis_cover(c(0, ymax))
  top <- ax$top
  y_model <- -0.10 * ymax   # schematic center line
  utr_half <- 0.018 * ymax  # UTR rectangle half-height
  cds_half <- 0.055 * ymax  # CDS rectangle half-height
  y_label <- -0.235 * ymax
  y_count <- ymax + 0.6 * (top - ymax)
  y_bot <- -0.31 * ymax

  # Boundary dashed lines stop at the y = 0 axis (geom_vline would cut
  # through the transcript schematic below)
  boundary_df <- tibble(x = region_bounds, xend = region_bounds, y = 0, yend = top)

  # Transcript schematic: three rectangles centered on y_model; thin UTRs and
  # a thick CDS
  model_df <- tibble(
    xmin = c(0, region_bounds[1], region_bounds[2]),
    xmax = c(region_bounds[1], region_bounds[2], 1),
    ymin = c(y_model - utr_half, y_model - cds_half, y_model - utr_half),
    ymax = c(y_model + utr_half, y_model + cds_half, y_model + utr_half)
  )

  region_label_df <- tibble(
    x = c(0.075, 0.5, 0.925),  # midpoints of the three segments
    y = y_label,
    label = c("5'UTR", "CDS", "3'UTR")
  )
  count_label_df <- curves |>
    transmute(Tissue, x = 0.5, y = y_count, label = sprintf("%d sites", n_sites))
  # Single y axis: the theme's axis.line.y would redraw the vertical line in
  # every panel, so draw it in the first panel only. Tissue must be a leveled
  # factor matching the main data -- a character vector would let
  # combine_vars reorder the panels alphabetically (Flower/Leaf/Root) and put
  # the axis line in the wrong panel
  axis_line_df <- tibble(Tissue = factor(tissue_levels[1], levels = tissue_levels),
                         x = 0, xend = 0, y = 0, yend = top)

  p <- ggplot() +
    geom_area(data = curve_df, aes(x, y, fill = Tissue)) +
    scale_fill_manual(values = tissue_colors, guide = "none") +
    geom_segment(data = boundary_df, aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "dashed", linewidth = 0.45) +
    geom_segment(data = axis_line_df, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = SIZE_TOKENS$base_linewidth) +
    annotate("segment", x = 0, xend = 1.04, y = 0, yend = 0,
             linewidth = SIZE_TOKENS$base_linewidth,
             arrow = arrow(length = unit(0.16, "cm"), type = "closed")) +
    geom_rect(data = model_df,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "black") +
    geom_text(data = region_label_df, aes(x, y, label = label), size = 16 / .pt) +
    geom_text(data = count_label_df, aes(x, y, label = label), size = 18 / .pt) +
    scale_x_continuous(limits = c(0, 1.08), breaks = NULL, expand = expansion(0, 0)) +
    scale_y_continuous(limits = c(y_bot, top), breaks = ax$breaks,
                       expand = expansion(mult = c(0.02, 0.05))) +
    facet_wrap(~Tissue,
               labeller = as_labeller(c(Root = "Roots", Leaf = "Leaves", Flower = "Flowers"))) +
    labs(x = NULL, y = "Density") +
    coord_cartesian(clip = "off") +
    theme(
      strip.background = element_blank(),
      panel.border = element_blank(),
      panel.spacing.x = grid::unit(20, "points"),
      axis.line = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = margin(10, 14, 8, 10)
    )

  outputs <- list(
    curve_tsv = file.path(out_dir, paste0(output_stem, "_curve.tsv")),
    plot_pdf = file.path(out_dir, paste0(output_stem, ".pdf")),
    plot_svg = file.path(out_dir, paste0(output_stem, ".svg")),
    plot_png = file.path(out_dir, paste0(output_stem, ".png"))
  )

  if (save_outputs) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(curve_df, outputs$curve_tsv)
    save_plot(p, file.path(out_dir, output_stem), width = 10.5, height = 4.8)
    message("Saved density curves to: ", outputs$curve_tsv)
    message("Saved plot files to: ", out_dir)
  }

  list(plot = p, curve_table = curve_df, output_paths = outputs, out_dir = out_dir)
}

if (sys.nframe() == 0) {
  invisible(run_metagene_density())
}
