#!/usr/bin/env Rscript
# plot.R — trua3 Root：23 个跨组织 dependent 位点的逐重复深度热图
#
# 展示目的：这批位点在 leaf/Flower 双株系都能判上 dependent，Root 里修饰明显下降
# 却全部卡在 q≥0.05 —— 根因是 WT Root 文库在 tRNA 位点上的覆盖度只有突变株的
# ~1/45，methylKit 加权 GLM 的过离散校正吃掉全部功效。本图按位点 × 样本
# （3 基因型 × 3 Treat 重复）铺开深度，格内标原始值，底部附每个样本在全部
# 439 个候选位点上的中位深度（说明失衡是文库级的，不限于这 23 个位点）。
#
# 输入（跑 data_prep.sh 生成）：
#   trua3_root_union23_replicate_depth.tsv  长表 site/fix_pos/genotype/rep/depth
#   trua3_root_candidate_median_depth.tsv   各样本在 439 候选上的中位深度

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
shared_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "..", "shared"), mustWork = TRUE)
source(file.path(shared_dir, "plot_style.R"))
DPI <- 600

# ComplexHeatmap 无设备测宽时的已知无害警告（与 ../dependent_trna_heatmap_common.R 同款）
muffle_arial_postscript_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("PostScript", msg) && grepl("Arial", msg)) invokeRestart("muffleWarning")
  })
}

text_gp <- function(size, face = "bold") {
  grid::gpar(fontsize = size, fontface = face, fontfamily = "Arial", col = "black")
}

long <- read.table(
  file.path(script_dir, "trua3_root_union23_replicate_depth.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)
median_df <- read.table(
  file.path(script_dir, "trua3_root_candidate_median_depth.tsv"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE
)

genotype_levels <- c("WT", "trua3_1", "trua3_2")
rep_levels <- c("rep1", "rep2", "rep3")
sample_levels <- as.vector(t(outer(genotype_levels, rep_levels, paste, sep = ".")))

mat_long <- long |>
  transmute(
    site = factor(
      paste0(sub("_[0-9]+$", "", sub("^tRNA-Nu-", "", site)), ":", sub(".*_", "", site)),
      levels = paste0(
        sub("_[0-9]+$", "", sub("^tRNA-Nu-", "", unique(long$site))),
        ":", sub(".*_", "", unique(long$site))
      )
    ),
    sample = factor(paste0(genotype, ".", rep), levels = sample_levels),
    depth
  )
mat <- mat_long |>
  tidyr::pivot_wider(names_from = sample, values_from = depth) |>
  tibble::column_to_rownames("site") |>
  as.matrix()
stopifnot(dim(mat) == c(23, 9))

log_mat <- log10(mat)
col_fun <- circlize::colorRamp2(
  c(1, 2, 3, 4, 5),
  c("#FFFFFF", "#FEE0D2", "#FC9272", "#DE2D26", "#A50F15")
)

# 格内标原始深度；底色深时改白字
label_mat <- matrix(sprintf("%.0f", mat), nrow = nrow(mat))
text_col_mat <- matrix(ifelse(log_mat > 3.3, "white", "black"), nrow = nrow(mat))

# 底部注释：439 候选上的中位深度（失衡是文库级）
fmt_depth <- function(v) ifelse(v >= 1000, paste0(round(v / 1000, 1), "k"), sprintf("%.0f", v))
median_vec <- median_df |>
  mutate(sample = paste0(sample, ".", rep)) |>
  arrange(match(sample, sample_levels)) |>
  pull(median_depth)
stopifnot(length(median_vec) == 9)

bottom_anno <- ComplexHeatmap::HeatmapAnnotation(
  Median = ComplexHeatmap::anno_text(
    fmt_depth(median_vec),
    gp = text_gp(8.5),
    rot = 0,
    show_name = TRUE
  ),
  annotation_name_side = "left",
  annotation_name_gp = text_gp(9),
  annotation_name_rot = 0
)

ht <- ComplexHeatmap::Heatmap(
  log_mat,
  name = "Depth",
  col = col_fun,
  bottom_annotation = bottom_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_order = seq_len(nrow(mat)),
  column_order = seq_len(ncol(mat)),
  row_names_side = "left",
  row_names_gp = text_gp(10.5),
  row_title = NULL,
  column_labels = rep(c("Rep 1", "Rep 2", "Rep 3"), 3),
  column_names_gp = text_gp(10),
  column_names_rot = 0,
  column_split = factor(rep(genotype_levels, each = 3), levels = genotype_levels),
  column_title = genotype_levels,
  column_title_gp = text_gp(14),
  column_gap = grid::unit(3, "mm"),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(
      label_mat[i, j], x, y,
      gp = grid::gpar(
        fontsize = 8.5, fontface = "bold", fontfamily = "Arial",
        col = text_col_mat[i, j]
      )
    )
  },
  heatmap_legend_param = list(
    title = "Depth",
    at = log10(c(30, 100, 300, 1000, 3000, 10000)),
    labels = c("30", "100", "300", "1k", "3k", "10k"),
    title_gp = text_gp(10.5),
    labels_gp = text_gp(10),
    grid_height = grid::unit(4.5, "mm"),
    grid_width = grid::unit(5, "mm"),
    legend_direction = "vertical"
  )
)

# 画布与字号（宽 8 in：行名/列名 10-14pt ≈ 比值 1.7-2.4，格内 8.5pt）
FIG_WIDTH <- 8
FIG_HEIGHT <- 7

muffle_arial_postscript_warnings(
  save_plot(
    plot = function() {
      grid::grid.newpage()
      ComplexHeatmap::draw(
        ht,
        heatmap_legend_side = "right",
        padding = grid::unit(c(4, 4, 4, 8), "mm")
      )
    },
    filename = file.path(script_dir, "trua3_root_union23_replicate_depth_heatmap"),
    width = FIG_WIDTH,
    height = FIG_HEIGHT,
    formats = c("pdf", "png", "svg"),
    dpi = DPI
  )
)
