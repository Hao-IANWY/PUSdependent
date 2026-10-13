#!/usr/bin/env Rscript

# mRNA metagene 位点密度图(metaPlotR 版式):
# 三个组织(Root/Leaf/Flower)的 WT 阳性位点投影到 5'UTR-CDS-3'UTR 归一化
# metagene 坐标(由 build_metagene_positions.sh 生成, 各区域固定占比 15%/70%/15%)
# 后画核密度曲线; 虚线为区域交界, 面板下方为转录本模式图(UTR 细矩形 + CDS 粗矩形,
# 两者垂直居中对齐)与区域标签, 虚线内侧顶部注释位点数, x 轴右端箭头。
# 三面板共用一根 y 轴(仅最左面板出刻度); 分面条带保留文字(复数命名
# Roots/Leaves/Flowers)但无边框底色, 不设图例。
#
# 输入: meta/mrna_metagene_sites.tsv (build_metagene_positions.sh 的输出)
# 输出: meta/mrna_metagene_density.{pdf,png,svg}
#       meta/mrna_metagene_density_curve.tsv (绘图的密度曲线数据)

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
  # meta 比 Figures 深一层: meta -> mrna_region_distribution -> wildtype_landscape -> Figures -> repo
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

  # 区域交界虚线位置, 与 build_metagene_positions.sh 的 U5_FRAC/CDS_FRAC 一致
  region_bounds <- c(0.15, 0.85)

  # 每组织核密度曲线(0-1 内截断不外溢); 单一 y 轴, 刻度按全局峰值罩顶,
  # 模式图/注释的 y 定位也统一按全局峰值的相对高度(三面板同坐标系)
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
  y_model <- -0.10 * ymax   # 模式图中心线
  utr_half <- 0.018 * ymax  # UTR 矩形半高
  cds_half <- 0.055 * ymax  # CDS 矩形半高
  y_label <- -0.235 * ymax
  y_count <- ymax + 0.6 * (top - ymax)
  y_bot <- -0.31 * ymax

  # 交界虚线只画到 y=0 轴线为止(geom_vline 会穿透下方基因模式图区域)
  boundary_df <- tibble(x = region_bounds, xend = region_bounds, y = 0, yend = top)

  # 转录本模式图: 三段均为以 y_model 为中心的矩形, UTR 细、CDS 粗
  model_df <- tibble(
    xmin = c(0, region_bounds[1], region_bounds[2]),
    xmax = c(region_bounds[1], region_bounds[2], 1),
    ymin = c(y_model - utr_half, y_model - cds_half, y_model - utr_half),
    ymax = c(y_model + utr_half, y_model + cds_half, y_model + utr_half)
  )

  region_label_df <- tibble(
    x = c(0.075, 0.5, 0.925),  # 三段区域的中点
    y = y_label,
    label = c("5'UTR", "CDS", "3'UTR")
  )
  count_label_df <- curves |>
    transmute(Tissue, x = 0.5, y = y_count, label = sprintf("%d sites", n_sites))
  # 单一 y 轴: theme 的 axis.line.y 会在每个面板重复画竖线, 改为只在首面板画。
  # Tissue 必须与主数据同为带 levels 的 factor —— 字符型会被 combine_vars 按字母序
  # 重排面板(Flower/Leaf/Root)且竖线落到错误的面板
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
