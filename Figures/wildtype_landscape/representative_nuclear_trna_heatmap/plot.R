#!/usr/bin/env Rscript
# plot.R — 代表性核 tRNA 的 positive-site ratio 沿 canonical tRNA 结构布局的热图。
#
# 行 = 46 个核 tRNA 家族，代表成员取自 tRNA/count/ncRNA/family_top/family_top_tRNA.tsv
#      （每个 anticodon 家族里 mean_CPM 最高的 tRNA）。行序：先按氨基酸字母序，
#      同氨基酸再按反密码子字母序；eMet/iMet 同属 Met，按 Met 归位（eMet 在前）。
# 列 = canonical 0–76 位（按 17 个区段分块）。
# 颜色 = 该家族在该 canonical 位上所有拷贝位点的 depth 加权平均 ratio。
# Sprinzl 坐标由 data.sh 直接取自 WT_*_sites.positive.Sprinzl.tsv，本脚本只做
# canonical 位号 -> 展示列的布局查表。
# 输入由 data.sh 生成（representative_nuclear_trna_mapped_sites_<tissue>.tsv）。
#
# 实现照搬 archive/20260902/paper/WT_sites_overview/representative_nuclear_trna_heatmap/
# plot_representative_nuclear_trna_heatmap.R,样式改走 Figures/shared。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
source(file.path(script_dir, "..", "..", "shared", "plot_style.R"))
source(file.path(script_dir, "..", "..", "shared", "plot_tokens.R"))

export_dpi <- 600
dimname_padding <- grid::unit(30, "mm")
heatmap_draw_padding <- grid::unit(c(4, 14, 4, 54), "mm")
ratio_palette <- c(low = "#FFF6EB", high = "#C44E00")

# 文字统一走 shared 的无 showtext 路线:gpar 里显式指定 Arial
text_gp <- function(size, color = "black") {
  grid::gpar(fontsize = size, fontface = "bold", fontfamily = "Arial", col = color)
}

# ============================================================================
# canonical 结构布局:区段顺序 / 固定长度 / 底部区段名 / 区段色
# ============================================================================
region_levels <- c(
  "p1", "Accstem_5p", "p8-9", "Dstem_5p", "Dloop", "Dstem_3p",
  "p26", "Acstem_5p", "Acloop", "Acstem_3p", "Vregion",
  "Tstem_5p", "Tloop", "Tstem_3p", "Accstem_3p", "p73", "CCA"
)

layout_config <- tibble::tibble(
  region = factor(region_levels, levels = region_levels),
  fixed_length = c(1L, 7L, 2L, 4L, 8L, 4L, 1L, 5L, 7L, 5L, 5L, 5L, 7L, 5L, 7L, 1L, 3L),
  broad_label = c(
    "", "Acceptor\nStem", "", "D Stem", "D Loop", "D Stem",
    "", "Anticodon\nStem", "Anticodon\nLoop", "Anticodon\nStem", "Variable",
    "T Stem", "T Loop", "T Stem", "Acceptor\nStem", "", "CCA"
  ),
  color_hex = c(
    "#FFFFFF", "#C93855", "#FFFFFF", "#EE7621", "#FFD700", "#EE7621",
    "#FFFFFF", "#F0E68C", "#FFFACD", "#F0E68C", "#90EE90",
    "#20B2AA", "#008B8B", "#20B2AA", "#C93855", "#FFFFFF", "#6A5ACD"
  )
)

column_metadata <- layout_config |>
  mutate(rel_pos = purrr::map(fixed_length, seq_len)) |>
  tidyr::unnest(rel_pos) |>
  mutate(
    matrix_col_id = paste(region, rel_pos, sep = "_"),
    display_number = 0:(dplyr::n() - 1L)
  )

matrix_col_ids <- column_metadata$matrix_col_id
max_canonical_pos <- max(column_metadata$display_number)
region_color_map <- setNames(layout_config$color_hex, as.character(layout_config$region))

# ============================================================================
# 读入映射结果,折叠成 代表 x canonical 位 的矩阵
# ============================================================================
read_tissue_matrix <- function(mapped_sites_path, ordered_rep_ids) {
  sites <- readr::read_tsv(mapped_sites_path, show_col_types = FALSE)

  required_cols <- c("Rep_ID", "canonical_pos", "ratio", "depth")
  missing_cols <- setdiff(required_cols, names(sites))
  if (length(missing_cols) > 0) {
    stop(basename(mapped_sites_path), " is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  unknown_reps <- setdiff(unique(sites$Rep_ID), ordered_rep_ids)
  if (length(unknown_reps) > 0) {
    stop("映射结果里出现了 family_top 之外的家族: ", paste(unknown_reps, collapse = ", "))
  }

  # canonical 位号 -> 展示列：纯布局查表（display_number 为 0..76）
  col_id_by_canonical <- column_metadata$matrix_col_id
  sites <- sites |>
    mutate(canonical_pos = as.integer(canonical_pos)) |>
    filter(
      !is.na(canonical_pos),
      canonical_pos >= 0L,
      canonical_pos <= max_canonical_pos
    ) |>
    mutate(matrix_col_id = col_id_by_canonical[canonical_pos + 1L])

  collapsed <- sites |>
    group_by(Rep_ID, matrix_col_id) |>
    summarise(
      ratio = weighted.mean(ratio, w = pmax(depth, 0), na.rm = TRUE),
      total_depth = sum(depth, na.rm = TRUE),
      n_sites = dplyr::n(),
      .groups = "drop"
    )

  # 一个 canonical 位可能收到多个原始位点(区段缩短时),已在上面按 depth 加权合并
  compressed <- collapsed |> filter(n_sites > 1)
  if (nrow(compressed) > 0) {
    message("collapsed ", nrow(compressed), " canonical cells from multiple source sites")
  }

  mat <- matrix(
    NA_real_,
    nrow = length(ordered_rep_ids),
    ncol = length(matrix_col_ids),
    dimnames = list(ordered_rep_ids, matrix_col_ids)
  )
  if (nrow(collapsed) > 0) {
    mat[cbind(match(collapsed$Rep_ID, ordered_rep_ids), match(collapsed$matrix_col_id, matrix_col_ids))] <-
      collapsed$ratio
  }
  mat
}

# ============================================================================
# 绘图
# ============================================================================
build_bottom_annotation <- function() {
  ComplexHeatmap::HeatmapAnnotation(
    Structure_Colors = ComplexHeatmap::anno_simple(
      x = as.character(column_metadata$region),
      col = region_color_map,
      height = grid::unit(3, "mm"),
      border = TRUE,
      gp = grid::gpar(col = "black", lwd = SIZE_TOKENS$base_linewidth)
    ),
    Numbers = ComplexHeatmap::anno_text(
      x = column_metadata$display_number,
      rot = 90,
      just = "right",
      gp = text_gp(SIZE_TOKENS$annotation_text_size)
    ),
    Regions = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = "white", col = NA),
      labels = layout_config$broad_label,
      labels_gp = text_gp(SIZE_TOKENS$annotation_text_size),
      labels_rot = 0
    ),
    gap = grid::unit(1, "mm"),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
}

build_ratio_legend <- function() {
  ComplexHeatmap::Legend(
    col_fun = circlize::colorRamp2(c(0, 1), c(ratio_palette[["low"]], ratio_palette[["high"]])),
    title = "Ratio",
    at = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    direction = "vertical",
    border = "black",
    legend_gp = grid::gpar(col = "black", lwd = SIZE_TOKENS$base_linewidth),
    background = "white",
    title_gp = text_gp(SIZE_TOKENS$legend_title_size),
    labels_gp = text_gp(SIZE_TOKENS$annotation_text_size),
    grid_width = grid::unit(7, "mm"),
    legend_height = grid::unit(54, "mm"),
    tick_length = grid::unit(3, "mm")
  )
}

build_tissue_heatmap <- function(mat) {
  old_dimname_padding <- ComplexHeatmap::ht_opt("DIMNAME_PADDING")
  on.exit(ComplexHeatmap::ht_opt(DIMNAME_PADDING = old_dimname_padding), add = TRUE)
  ComplexHeatmap::ht_opt(DIMNAME_PADDING = dimname_padding)

  ComplexHeatmap::Heatmap(
    mat,
    name = "Ratio",
    col = circlize::colorRamp2(c(0, 1), c(ratio_palette[["low"]], ratio_palette[["high"]])),
    na_col = "#FCFCFC",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = factor(column_metadata$region, levels = layout_config$region),
    column_gap = grid::unit(0, "mm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = SIZE_TOKENS$base_linewidth),
    bottom_annotation = build_bottom_annotation(),
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_names_gp = text_gp(SIZE_TOKENS$annotation_text_size),
    row_names_max_width = grid::unit(8, "cm"),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
      title_gp = text_gp(SIZE_TOKENS$legend_title_size),
      labels_gp = text_gp(SIZE_TOKENS$annotation_text_size)
    )
  )
}

# 热图本体 + 右侧带边框的竖直图例(图例由脚本手工放置,与参考图一致)
draw_tissue_heatmap <- function(ht) {
  legend_obj <- build_ratio_legend()
  legend_width <- ComplexHeatmap:::width(legend_obj) + grid::unit(6, "mm")
  legend_height <- ComplexHeatmap:::height(legend_obj) + grid::unit(6, "mm")

  ComplexHeatmap::draw(
    ht,
    merge_legends = FALSE,
    show_heatmap_legend = FALSE,
    show_annotation_legend = FALSE,
    padding = heatmap_draw_padding,
    newpage = FALSE
  )

  grid::upViewport(0)
  grid::pushViewport(
    grid::viewport(
      x = grid::unit(1, "npc") - heatmap_draw_padding[4] / 2,
      y = grid::unit(0.5, "npc"),
      width = legend_width,
      height = legend_height,
      just = c("center", "center"),
      clip = "off"
    )
  )
  grid::grid.rect(gp = grid::gpar(fill = NA, col = "black", lwd = SIZE_TOKENS$base_linewidth))
  ComplexHeatmap::draw(
    legend_obj,
    x = grid::unit(2, "mm"),
    y = grid::unit(0.5, "npc"),
    just = c("left", "center")
  )
  grid::popViewport()
}

# ============================================================================
# 主流程
# ============================================================================
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
family_top_path <- file.path(repo_root, "tRNA/count/ncRNA/family_top/family_top_tRNA.tsv")
if (!file.exists(family_top_path)) {
  stop("Missing family_top table: ", family_top_path)
}
family_top <- readr::read_tsv(family_top_path, show_col_types = FALSE)

# 行 = 核 tRNA 家族；代表成员（mean_CPM 最高）记录在 Representative_ID 列
# 行序（老板 2026-10-25 定）：先按氨基酸字母序，同氨基酸再按反密码子字母序；
# eMet/iMet 是小写开头，若按整名排序会沉底/散开，故归一成 Met 参与排序
ordered_rep_ids <- family_top |>
  filter(str_starts(feature_id, "^tRNA-Nu-")) |>
  mutate(
    aa = str_replace(str_remove(str_remove(family, "^Nu-"), "-[ACGT]{3}$"),
                     "^[ei]Met$", "Met"),
    anticodon = str_extract(family, "[ACGT]{3}$")
  ) |>
  arrange(aa, anticodon, family) |>
  pull(family)
if (length(ordered_rep_ids) == 0) {
  stop("family_top 表里没有核 tRNA 家族: ", family_top_path)
}

tissues <- tibble::tribble(
  ~Tissue, ~prefix,
  "Root", "root",
  "Leaf", "leaf",
  "Flower", "flower"
)

for (i in seq_len(nrow(tissues))) {
  tissue <- tissues$Tissue[i]
  prefix <- tissues$prefix[i]
  mapped_path <- file.path(script_dir, paste0("representative_nuclear_trna_mapped_sites_", prefix, ".tsv"))
  if (!file.exists(mapped_path)) {
    stop("Missing mapped sites for ", tissue, ": ", mapped_path, "\n先跑 data.sh 生成。")
  }

  mat <- read_tissue_matrix(mapped_path, ordered_rep_ids)
  n_cells <- sum(!is.na(mat))
  message(tissue, ": ", length(ordered_rep_ids), " 行 x ", length(matrix_col_ids),
          " 个 canonical 位, 有值单元格 ", n_cells)

  # 行数固定 46,高度公式与参考图一致
  height_in <- max(6.5, 0.18 * length(ordered_rep_ids) + 2.8)

  # ComplexHeatmap 在"无设备"状态下测量文本宽度时,R 会临时打开默认 pdf 设备,而该设备的
  # PostScript 字体库里没有 Arial,于是每个被测量的文本元素抛一条警告。测量值不受影响
  # (实测与 cairo 下一致),且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。
  # 这里只屏蔽这一条已知无害的警告,其他警告照常抛出。
  withCallingHandlers({
    ht <- build_tissue_heatmap(mat)
    save_plot(
      plot = function() {
        grid::grid.newpage()
        draw_tissue_heatmap(ht)
      },
      filename = file.path(script_dir, paste0("representative_nuclear_trna_heatmap_", prefix)),
      width = 18, height = height_in, formats = c("pdf", "png", "svg"), dpi = export_dpi
    )
  }, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("PostScript", msg) && grepl("Arial", msg)) invokeRestart("muffleWarning")
  })
}

message("[plot] Saved representative nuclear tRNA heatmaps to: ", script_dir)
