#!/usr/bin/env Rscript
# ============================================================================
# plot.R — 「27℃ 温度上调位点」vs「27℃ trua1b dependent 位点」overlap Venn
#          （mRNA / rRNA / tRNA，各一张）
#
# 输入：<rna>_27up_vs_trua1b_dependent_overlap.tsv（data_prep.sh 生成）
# 输出：<rna>_27up_vs_trua1b_dependent_venn.pdf / .png / .svg
#   左圈 = 27℃ 上调位点，右圈 = 27℃ trua1b dependent 位点，中间是交集计数。
#
# 移植自论文版 archive/20260902/paper/TRUA1B_thermomorphogenesis/
# temperature_responsive_dependent_overlap/plot.R 的 "27only" 2-set 变体
# （原图只有 mRNA 一张），这里拆成三个 RNA 类型，样式改走 Figures/shared。
#
# 注意：本分析在新版数据里三组交集全为 0（见 README）。用的是手绘 ggplot 圆
# 而非 VennDiagram——手绘的好处正是空集合也能照常出图，不会因交集为 0 报错。
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

register_font()

rna_types <- c("mRNA", "rRNA", "tRNA")

set_up <- "27°C up-regulated"
set_dependent <- "TRUA1B-dependent at 27°C"
# 二色与 archive 的 2-set 变体一致（alpha 0x66），属本图专用，不入 plot_tokens
venn_colors <- c("#E64B3566", "#3C548866")
names(venn_colors) <- c(set_up, set_dependent)

venn_text_size <- SIZE_TOKENS$axis_title_size / ggplot2::.pt

# ggplot 画圆：VennDiagram 那类包不支持空集合，手绘圆不受集合大小/交集影响
circle_df <- function(cx, cy, r, label, n = 360) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(
    x = cx + r * cos(theta),
    y = cy + r * sin(theta),
    set = label
  )
}

# 图内不留坐标轴/网格，只保留 Arial bold 的文本
venn_theme <- function() {
  theme_void() +
    theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      plot.margin = margin(10, 10, 10, 10)
    )
}

build_venn <- function(up_only, dependent_only, both) {
  circles <- bind_rows(
    circle_df(-0.75, 0, 1.5, set_up),
    circle_df(0.75, 0, 1.5, set_dependent)
  )

  counts <- tibble(
    x = c(-1.5, 1.5, 0.0),
    y = c(0.0, 0.0, 0.0),
    label = c(up_only, dependent_only, both)
  )

  labels <- tibble(
    x = c(-1.5, 1.5),
    y = c(2.0, 2.0),
    label = c(
      "bold('27'*degree*plain(C)~'up-regulated')",
      "bold('TRUA1B-dependent at '*27*degree*plain(C))"
    )
  )

  ggplot() +
    geom_polygon(
      data = circles,
      aes(x = .data$x, y = .data$y, group = .data$set, fill = .data$set),
      color = "black",
      linewidth = 0.8
    ) +
    geom_text(
      data = counts,
      aes(x = .data$x, y = .data$y, label = .data$label),
      size = venn_text_size,
      family = "Arial",
      fontface = "bold",
      colour = "black"
    ) +
    geom_text(
      data = labels,
      aes(x = .data$x, y = .data$y, label = .data$label),
      size = venn_text_size,
      family = "Arial",
      fontface = "bold",
      parse = TRUE
    ) +
    scale_fill_manual(values = venn_colors, guide = "none") +
    coord_equal(xlim = c(-2.8, 2.8), ylim = c(-2.0, 2.5), clip = "off") +
    venn_theme()
}

for (rna_type in rna_types) {
  stem <- paste0(rna_type, "_27up_vs_trua1b_dependent_venn")
  input_file <- file.path(script_dir, paste0(rna_type, "_27up_vs_trua1b_dependent_overlap.tsv"))
  if (!file.exists(input_file)) {
    stop("Input table not found: ", input_file, "\nRun data_prep.sh first.")
  }

  overlap_df <- read.table(
    input_file,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(overlap_df) == 0) {
    stop("No rows found in overlap table: ", input_file)
  }

  combo_counts <- table(factor(overlap_df$combination,
                               levels = c("up_only", "dependent_only", "both")))
  up_only <- as.integer(combo_counts[["up_only"]])
  dependent_only <- as.integer(combo_counts[["dependent_only"]])
  both <- as.integer(combo_counts[["both"]])

  if (both == 0) {
    message("[plot] ", rna_type, ": overlap is 0 (", up_only, " up-only vs ",
            dependent_only, " dependent-only); the intersection cell shows 0.")
  }

  save_plot(
    build_venn(up_only, dependent_only, both),
    file.path(script_dir, stem),
    width = 9,
    height = 7
  )
  message(sprintf("[plot] %s: up_only=%d, dependent_only=%d, both=%d",
                  rna_type, up_only, dependent_only, both))
}

message("[plot] Done: 3 temperature-responsive vs trua1b-dependent overlap Venns (pdf/png/svg)")
