#!/usr/bin/env Rscript
# ============================================================================
# plot.R — call_dependent_methylKit.R 统计学检验流程示意图
#
# 纯概念流程图（无数据来源），五个步骤自上而下：
#   1 输入：一个候选位点的 6 个文库 (depth, ratio) → 计数
#   2 模型：加权二项 logistic 回归 logit(π) = β0 + β1·x
#   3 检验：过离散 φ 校正（φ>1 走 F，φ≤1 走卡方）→ 双侧 p
#   4 多重校正：在所有候选位点上做 BH 校正 → q 值
#   5 判定：方向为正（meth.diff > 0）、mock depth ≥ 10、q < 0.05 且 Δ ≥ 0.05，
#     两株系取交集
#
# 输出：dependent_methylKit_test_schematic.pdf / .png / .svg
# 原理文档：同目录 README.md
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))

register_font()

# ---- 配色（沿用 PAL_METHOD 语义：蓝 = WT/dependent 体系，红 = 突变/methylKit 强调） ----
COL_WT   <- "#3C5488"
COL_MUT  <- "#E64B35"
COL_GREY <- "grey35"

# ---- 绘图助手（显式 family = "Arial"：geom_text 缺省会落到 DejaVu） ----
box <- function(x0, x1, y0, y1, label, border = "black", fill = "white",
                size = 3.1, dashed = FALSE, fontface = "plain",
                colour = "black", lineheight = 0.95) {
  list(
    annotate("rect", xmin = x0, xmax = x1, ymin = y0, ymax = y1,
             fill = fill, colour = border, linewidth = 0.7,
             linetype = if (dashed) "dashed" else "solid"),
    annotate("text", x = (x0 + x1) / 2, y = (y0 + y1) / 2, label = label,
             size = size, family = "Arial", fontface = fontface,
             colour = colour, lineheight = lineheight)
  )
}

tbox <- function(x, y, label, size = 3.1, fontface = "plain", colour = "black",
                 hjust = 0.5, lineheight = 0.95) {
  annotate("text", x = x, y = y, label = label, size = size, family = "Arial",
           fontface = fontface, colour = colour, hjust = hjust,
           lineheight = lineheight)
}

arrow_down <- function(x, y0, y1) {
  annotate("segment", x = x, xend = x, y = y0, yend = y1,
           arrow = arrow(length = unit(0.09, "in")), linewidth = 0.7,
           colour = "grey25")
}

arrow_to <- function(x0, y0, x1, y1) {
  annotate("segment", x = x0, y = y0, xend = x1, yend = y1,
           arrow = arrow(length = unit(0.09, "in")), linewidth = 0.7,
           colour = "grey25")
}

light <- function(col) paste0(col, "1A")

p <- ggplot() +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 102), expand = FALSE) +
  theme_void() +

  # ---- 标题 ----
  tbox(50, 99.5, "How call_dependent_methylKit.R tests one candidate site",
       size = 4.5, fontface = "bold") +

  # ---- 1 输入：6 个文库 ----
  tbox(2, 96.5, "1 | Input — one candidate site, six libraries",
       size = 3.8, fontface = "bold", hjust = 0) +
  tbox(27.5, 93.2, "WT  (x = 1)", size = 3.4, fontface = "bold", colour = COL_WT) +
  tbox(72.5, 93.2, "mutant  (x = 0)", size = 3.4, fontface = "bold", colour = COL_MUT) +
  box( 6, 19, 85.5, 91, "rep1\n(depth, ratio)", border = COL_WT,  fill = light(COL_WT),  size = 3.0) +
  box(21, 34, 85.5, 91, "rep2\n(depth, ratio)", border = COL_WT,  fill = light(COL_WT),  size = 3.0) +
  box(36, 49, 85.5, 91, "rep3\n(depth, ratio)", border = COL_WT,  fill = light(COL_WT),  size = 3.0) +
  box(51, 64, 85.5, 91, "rep1\n(depth, ratio)", border = COL_MUT, fill = light(COL_MUT), size = 3.0) +
  box(66, 79, 85.5, 91, "rep2\n(depth, ratio)", border = COL_MUT, fill = light(COL_MUT), size = 3.0) +
  box(81, 94, 85.5, 91, "rep3\n(depth, ratio)", border = COL_MUT, fill = light(COL_MUT), size = 3.0) +
  arrow_down(50, 85.5, 83.2) +

  # 计数还原
  box(20, 80, 76.5, 83,
      "each sample → counts:  modified  k = round(depth × ratio),  unmodified  depth − k\ndepth matters: 10/20 carries more weight than 1/2",
      size = 3.0) +
  arrow_down(50, 76.5, 74.7) +

  # 零覆盖处理（预处理，虚线）
  box(8, 92, 67.5, 74.5,
      "depth = 0 replicate  →  treated as MISSING (not \"modification lost\")\nexcluded from Mut_Ratio mean and from the test:  3v3 → 3v2\nif a whole group is wiped out  →  p = NA",
      border = COL_GREY, dashed = TRUE, size = 2.9) +
  arrow_down(50, 67.5, 65.7) +

  # ---- 2 模型 ----
  tbox(2, 64, "2 | Weighted binomial logistic regression (per site)",
       size = 3.8, fontface = "bold", hjust = 0) +
  box(12, 88, 56.5, 63,
      "logit(π) = β0 + β1 · x        x = 1 (WT) / 0 (mutant),   weight = depth\nβ1 > 0  means  WT > mutant   →   the dependent direction",
      size = 3.2) +
  arrow_down(50, 56.5, 54.7) +

  # ---- 3 过离散校正检验 ----
  tbox(2, 53, "3 | Overdispersion-corrected test",
       size = 3.8, fontface = "bold", hjust = 0) +
  box(30, 70, 48.5, 52, "overdispersion factor  φ  (MN)", size = 3.2) +
  arrow_to(44, 48.5, 28, 44.8) +
  arrow_to(56, 48.5, 72, 44.8) +
  box(6, 48, 37.5, 44.5,
      "φ > 1   →   F test\ndeviance/φ ~ F(1, n−2)\npenalizes noisy replicates",
      size = 2.9) +
  box(52, 94, 37.5, 44.5,
      "φ ≤ 1   →   χ² test\nordinary likelihood ratio\nreplicates consistent",
      size = 2.9) +
  arrow_to(27, 37.5, 46, 35) +
  arrow_to(73, 37.5, 54, 35) +
  box(30, 70, 31, 34.8,
      "two-sided p-value\n→  BH q-value across all candidate sites",
      size = 3.0) +
  arrow_down(50, 31, 27.2) +

  # ---- 4 方向 + FDR ----
  tbox(2, 28.4, "4 | Direction + FDR",
       size = 3.8, fontface = "bold", hjust = 0) +
  box(12, 88, 21.5, 27,
      "meth.diff > 0 (WT higher)   AND   q < 0.05   →   passes\nwrong direction or q ≥ 0.05   →   never passes",
      size = 3.1) +
  arrow_down(50, 21.5, 17.7) +

  # ---- 5 判定 ----
  tbox(2, 18.7, "5 | Decision gates — all four must hold",
       size = 3.8, fontface = "bold", hjust = 0) +
  list(
    annotate("rect", xmin = 12, xmax = 88, ymin = 9.5, ymax = 17.5,
             fill = "white", colour = "black", linewidth = 0.7),
    tbox(15, 15.9, "•  mock mean depth ≥ 10  —  site well measured in mutant mock",
         size = 3.0, hjust = 0),
    tbox(15, 14.1, "•  BH q < 0.05  —  significant loss in mutant (FDR)",
         size = 3.0, hjust = 0),
    tbox(15, 12.3, "•  Δ = WT_Ratio − Mut_Ratio ≥ 0.05  —  effect-size floor",
         size = 3.0, hjust = 0),
    tbox(15, 10.5, "•  both mutant lines pass all above  —  line 1 AND line 2",
         size = 3.0, hjust = 0)
  ) +
  arrow_down(50, 9.5, 7.7) +

  box(30, 70, 2.5, 7.5, "Dependent site",
      border = COL_MUT, fill = paste0(COL_MUT, "26"),
      size = 4.0, fontface = "bold")

save_plot(p, file.path(script_dir, "dependent_methylKit_test_schematic"),
          width = 9, height = 11)
