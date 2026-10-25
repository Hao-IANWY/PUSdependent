#!/usr/bin/env Rscript
# ============================================================================
# plot.R — WT 21/27 度温度响应 Ψ 位点散点图（mRNA / rRNA / tRNA 各一张）
#
# 输入：<rna>_temperature_responsive_scatter.tsv（data_prep.sh 生成）
#   每行一个候选位点（两个温度 fisher 阳性位点并集），带两个温度各自的修饰比
#   与 regulation（non-responsive / down-regulated / up-regulated）。
# 输出：<rna>_temperature_responsive_scatter.pdf / .png / .svg
#   X 轴 = WT 21 度修饰比例，Y 轴 = WT 27 度修饰比例，虚线为 y=x；
#   灰点 = 非温度响应位点，红点 = 27 度上调，蓝点 = 27 度下调。
#
# 移植自论文版 archive/20260902/paper/TRUA1B_thermomorphogenesis/
# temperature_responsive_sites/plot.R（原图只有 mRNA），这里拆成三张，
# 样式改走 Figures/shared。
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

rna_types <- c("mRNA", "rRNA", "tRNA")

ratio_breaks <- seq(0, 1, by = 0.2)
# 三个图共用同一套图例顺序与配色；drop = FALSE 保证某个类别空集时（tRNA 没有
# down-regulated 位点）图例仍然三行齐全，三张图可直接并排比较。
regulation_levels <- c("non-responsive", "down-regulated", "up-regulated")
regulation_colors <- c(
  "non-responsive" = "#B8B8B8CC",
  "down-regulated" = "#3C5488FF",
  "up-regulated"   = "#E64B35FF"
)

init_theme(base_size = SIZE_TOKENS$base_size, sizes = SIZE_TOKENS)

build_plot <- function(df) {
  # 灰点先画、响应位点后画，保证响应位点压在灰点之上
  ggplot(df, aes(x = .data$WT_21_Ratio, y = .data$WT_27_Ratio)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = SIZE_TOKENS$base_linewidth,
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_point(
      data = filter(df, .data$regulation == "non-responsive"),
      aes(fill = .data$regulation),
      shape = 21, size = 3.0, color = "black", stroke = 0.5, alpha = 0.88
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

message("[plot] Done: mRNA / rRNA / tRNA temperature-responsive scatter plots (pdf/png/svg)")
