#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggprism)
  library(sysfonts)
  library(showtext)
  library(scales)
  library(patchwork)
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

read_positive_sites <- function(file_path, tissue, rna_type) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "ratio")
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat <- dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      ratio_pct = ratio * 100
    )

  if (rna_type == "tRNA") {
    dat <- dat %>%
      filter(stringr::str_detect(chr, "^tRNA-"))
  }

  dat %>%
    filter(is.finite(ratio_pct), ratio_pct >= 0, ratio_pct <= 100) %>%
    transmute(
      tissue = factor(tissue, levels = c("Root", "Leaf", "Flower")),
      rna_type = factor(rna_type, levels = c("mRNA", "rRNA", "tRNA")),
      chr,
      pos,
      strand,
      site_id = paste(chr, pos, strand, sep = ":"),
      ratio_pct
    )
}

build_density_plot <- function(site_tbl, out_dir, output_stem) {
  panel_levels <- c(
    "mRNA-Root", "mRNA-Leaf", "mRNA-Flower",
    "rRNA-Root", "rRNA-Leaf", "rRNA-Flower",
    "tRNA-Root", "tRNA-Leaf", "tRNA-Flower"
  )

  site_tbl <- site_tbl %>%
    mutate(panel = factor(paste(as.character(rna_type), as.character(tissue), sep = "-"), levels = panel_levels))

  count_tbl <- site_tbl %>%
    count(tissue, rna_type, panel, name = "n_sites") %>%
    mutate(label = paste0("n = ", scales::comma(n_sites)))

  density_tbl <- site_tbl %>%
    group_by(tissue, rna_type, panel) %>%
    summarise(
      peak_density = {
        x <- ratio_pct[is.finite(ratio_pct)]
        if (length(x) < 2 || length(unique(x)) < 2) {
          0.02
        } else {
          max(stats::density(x, from = 0, to = 100, n = 512)$y, na.rm = TRUE)
        }
      },
      .groups = "drop"
    ) %>%
    left_join(count_tbl, by = c("tissue", "rna_type", "panel")) %>%
    mutate(
      tissue = factor(as.character(tissue), levels = c("Root", "Leaf", "Flower")),
      rna_type = factor(as.character(rna_type), levels = c("mRNA", "rRNA", "tRNA")),
      panel = factor(as.character(panel), levels = panel_levels)
    )

  fill_values <- c(
    "rRNA" = "#00A087CC",
    "tRNA" = "#3C5488CC",
    "mRNA" = "#E69F00CC"
  )

  density_plot <- ggplot(site_tbl, aes(x = ratio_pct, fill = rna_type)) +
    geom_density(
      color = "black",
      linewidth = SIZE_TOKENS$base_linewidth,
      alpha = 0.55,
      trim = TRUE
    ) +
    geom_text(
      data = density_tbl,
      aes(x = 50, y = peak_density * 0.95, label = label),
      inherit.aes = FALSE,
      family = "Arial",
      fontface = "bold",
      size = SIZE_TOKENS$annotation_text_size / ggplot2::.pt,
      hjust = 0.5,
      vjust = 1
    ) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    facet_wrap(~ panel, ncol = 3, scales = "free_y",
               # 内部 panel 键保持单数, 仅条带显示层改复数
               labeller = as_labeller(c(
                 "mRNA-Root" = "mRNA-Roots", "mRNA-Leaf" = "mRNA-Leaves", "mRNA-Flower" = "mRNA-Flowers",
                 "rRNA-Root" = "rRNA-Roots", "rRNA-Leaf" = "rRNA-Leaves", "rRNA-Flower" = "rRNA-Flowers",
                 "tRNA-Root" = "tRNA-Roots", "tRNA-Leaf" = "tRNA-Leaves", "tRNA-Flower" = "tRNA-Flowers"
               ))) +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, 20),
      expand = expansion(mult = c(0, 0)),
      name = "Modification ratio (%)"
    ) +
    scale_y_continuous(
      name = "Density",
      expand = expansion(mult = c(0, 0.06))
    ) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      panel.spacing.x = grid::unit(22, "pt"),
      panel.spacing.y = grid::unit(14, "pt"),
      strip.background = element_blank(),
      strip.placement = "outside",
      axis.line = element_line(color = "black", linewidth = SIZE_TOKENS$base_linewidth),
      plot.margin = margin(18, 24, 18, 24)
    )

  density_outputs <- list(
    sites_tsv = file.path(out_dir, paste0(output_stem, "_density_sites.tsv")),
    plot_pdf = file.path(out_dir, paste0(output_stem, "_density_raincloud.pdf")),
    plot_svg = file.path(out_dir, paste0(output_stem, "_density_raincloud.svg")),
    plot_png = file.path(out_dir, paste0(output_stem, "_density_raincloud.png"))
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(site_tbl, density_outputs$sites_tsv)
  save_plot(density_plot, file.path(out_dir, paste0(output_stem, "_density_raincloud")), width = 12.8, height = 8.8)

  list(plot = density_plot, output_paths = density_outputs, count_table = count_tbl)
}

build_overlap_table <- function(site_tbl, tissue_a, tissue_b, comparison_label) {
  a_tbl <- site_tbl %>%
    filter(as.character(tissue) == tissue_a) %>%
    select(rna_type, site_id, ratio_a = ratio_pct)

  b_tbl <- site_tbl %>%
    filter(as.character(tissue) == tissue_b) %>%
    select(rna_type, site_id, ratio_b = ratio_pct)

  inner_join(a_tbl, b_tbl, by = c("rna_type", "site_id")) %>%
    mutate(
      comparison = factor(comparison_label, levels = c("Root vs Leaf", "Root vs Flower", "Leaf vs Flower")),
      tissue_a = tissue_a,
      tissue_b = tissue_b
    )
}

build_correlation_plot <- function(site_tbl, out_dir, output_stem) {
  corr_site_tbl <- site_tbl %>%
    filter(rna_type %in% c("rRNA", "tRNA"))

  comparison_specs <- tribble(
    ~tissue_a, ~tissue_b, ~comparison,
    "Root", "Leaf", "Root vs Leaf",
    "Root", "Flower", "Root vs Flower",
    "Leaf", "Flower", "Leaf vs Flower"
  )

  overlap_tbl <- pmap_dfr(comparison_specs, build_overlap_table, site_tbl = corr_site_tbl)

  if (nrow(overlap_tbl) == 0) {
    stop("No overlap sites were found across the tissue comparisons.")
  }

  safe_r2 <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]
    if (length(x) < 2 || length(unique(x)) < 2 || length(unique(y)) < 2) {
      return(NA_real_)
    }
    stats::cor(x, y)^2
  }

  summary_tbl <- overlap_tbl %>%
    group_by(rna_type, comparison) %>%
    summarise(
      n_sites = n(),
      r2 = safe_r2(ratio_a, ratio_b),
      .groups = "drop"
    ) %>%
    mutate(
      rna_type = factor(as.character(rna_type), levels = c("rRNA", "tRNA")),
      comparison = factor(as.character(comparison), levels = c("Root vs Leaf", "Root vs Flower", "Leaf vs Flower")),
      label = sprintf("italic(r)^2 == %.3f", r2)
    )

  label_tbl <- summary_tbl %>%
    mutate(
      x = 10,
      y = 94
    )

  build_comparison_plot <- function(comp_name, rna_name, x_label, y_label) {
    plot_data <- overlap_tbl %>%
      filter(comparison == comp_name, rna_type == rna_name)

    plot_labels <- label_tbl %>%
      filter(comparison == comp_name, rna_type == rna_name)

    p <- ggplot(plot_data, aes(x = ratio_a, y = ratio_b)) +
      geom_point(
        shape = 21,
        size = 2.7,
        stroke = 0.5,
        color = "black",
        fill = "#D9D9D9",
        alpha = 0.95
      ) +
      geom_smooth(method = "lm", se = FALSE, linewidth = SIZE_TOKENS$base_linewidth, color = "#D95F02") +
      geom_text(
        data = plot_labels,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        family = "Arial",
        fontface = "bold",
        size = SIZE_TOKENS$annotation_text_size / ggplot2::.pt,
        hjust = 0,
        vjust = 1,
        parse = TRUE
      ) +
      scale_x_continuous(
        breaks = seq(0, 100, 25),
        expand = expansion(mult = c(0, 0))
      ) +
      scale_y_continuous(
        breaks = seq(0, 100, 25),
        expand = expansion(mult = c(0, 0))
      ) +
      coord_equal(xlim = c(0, 100), ylim = c(0, 100), clip = "off") +
      labs(
        x = x_label,
        y = y_label
      ) +
      theme(
        plot.margin = margin(30, 30, 32, 34),
        panel.border = element_blank(),
        axis.line = element_line(color = "black", linewidth = SIZE_TOKENS$base_linewidth)
      )

    p
  }

  corr_plot <- wrap_plots(
    # 第 1/2 参数是内部键(单数), 第 3/4 参数是显示在坐标轴上的组织名(复数)
    build_comparison_plot("Root vs Leaf", "rRNA", "Roots", "Leaves"),
    build_comparison_plot("Root vs Flower", "rRNA", "Roots", "Flowers"),
    build_comparison_plot("Leaf vs Flower", "rRNA", "Leaves", "Flowers"),
    build_comparison_plot("Root vs Leaf", "tRNA", "Roots", "Leaves"),
    build_comparison_plot("Root vs Flower", "tRNA", "Roots", "Flowers"),
    build_comparison_plot("Leaf vs Flower", "tRNA", "Leaves", "Flowers"),
    ncol = 3,
    byrow = TRUE
  ) &
    theme(
      legend.position = "none",
      plot.margin = margin(10, 10, 10, 10)
    )

  corr_outputs <- list(
    sites_tsv = file.path(out_dir, paste0(output_stem, "_correlation_sites.tsv")),
    summary_tsv = file.path(out_dir, paste0(output_stem, "_correlation_summary.tsv")),
    plot_pdf = file.path(out_dir, paste0(output_stem, "_correlation.pdf")),
    plot_svg = file.path(out_dir, paste0(output_stem, "_correlation.svg")),
    plot_png = file.path(out_dir, paste0(output_stem, "_correlation.png"))
  )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(overlap_tbl, corr_outputs$sites_tsv)
  readr::write_tsv(summary_tbl, corr_outputs$summary_tsv)
  save_plot(corr_plot, file.path(out_dir, paste0(output_stem, "_correlation")), width = 11.6, height = 7.8)

  list(plot = corr_plot, output_paths = corr_outputs, summary_table = summary_tbl)
}

run_analysis <- function(out_dir = NULL, output_stem = "rna_site_distribution_overlap") {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)
  out_dir <- if (is.null(out_dir)) script_dir else normalizePath(out_dir, mustWork = FALSE)

  input_files <- tribble(
    ~tissue, ~rna_type, ~file_path,
    "Root", "rRNA", file.path(repo_root, "polyARNA/STAT/fisher/rRNA/WT-Root_sites.positive.tsv"),
    "Leaf", "rRNA", file.path(repo_root, "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv"),
    "Flower", "rRNA", file.path(repo_root, "polyARNA/STAT/fisher/rRNA/WT-Flower_sites.positive.tsv"),
    "Root", "mRNA", file.path(repo_root, "polyARNA/STAT/fisher/mRNA/WT-Root_sites.positive.tsv"),
    "Leaf", "mRNA", file.path(repo_root, "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv"),
    "Flower", "mRNA", file.path(repo_root, "polyARNA/STAT/fisher/mRNA/WT-Flower_sites.positive.tsv"),
    "Root", "tRNA", file.path(repo_root, "tRNA/STAT/fisher/ncRNA/WT_Root_sites.positive.tsv"),
    "Leaf", "tRNA", file.path(repo_root, "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv"),
    "Flower", "tRNA", file.path(repo_root, "tRNA/STAT/fisher/ncRNA/WT_Flower_sites.positive.tsv")
  )

  missing_files <- input_files$file_path[!file.exists(input_files$file_path)]
  if (length(missing_files) > 0) {
    stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
  }

  site_tbl <- pmap_dfr(input_files, read_positive_sites)

  density_result <- build_density_plot(site_tbl, out_dir = out_dir, output_stem = output_stem)
  corr_result <- build_correlation_plot(site_tbl, out_dir = out_dir, output_stem = output_stem)

  message("Saved density plot to: ", density_result$output_paths$plot_pdf)
  message("Saved correlation plot to: ", corr_result$output_paths$plot_pdf)

  invisible(list(
    site_table = site_tbl,
    density = density_result,
    correlation = corr_result,
    out_dir = out_dir,
    output_stem = output_stem
  ))
}

if (sys.nframe() == 0) {
  run_analysis()
}
