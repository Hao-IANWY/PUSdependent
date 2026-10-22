#!/usr/bin/env Rscript
# dependent_trna_heatmap_common.R — 不同基因型 tRNA dependent 位点 delta-ratio 热图公共实现
#
# 移植自 archive/20260902/paper/tRNA_dependent/dependent_trna_heatmap_common.R，
# 样式改走 Figures/shared（plot_style.R + plot_tokens.R），导出改走 save_plot()。
#
# 两张图，共用同一套 canonical 列布局：
#   run_tissue_split_heatmap()   行 = 组织 × 代表 tRNA，高度随行数自适应（宽固定 18 in）
#   run_union_dependent_heatmap() 行 = 代表 tRNA，值 = 跨组织最大 Delta Ration
#
# 列布局 = canonical 0-76 位（p1 ... CCA 共 17 个区段、77 列）：每个区段按固定长度
# 重标定，底部有区段色带、canonical 位号、区段名。组织拆分图与 union 图均 Cy/Mt/Pt
# 三分室进图（2026-10-12 老板拍板；Pt tRNA 环长分布与核/线粒体不同，投影到同一布局
# 有扭曲，按既定口径接受）。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(purrr)
  library(stringr)
  library(tibble)
  library(grid)
})

export_dpi <- 600
dimname_padding <- grid::unit(24, "mm")
heatmap_draw_padding <- grid::unit(c(4, 16, 4, 54), "mm")
default_value_col <- "Delta_Ratio(WT-Mut)"
default_value_title <- "Delta Ratio"
tissue_levels <- c("Root", "leaf", "Flower")
tissue_display_map <- c(Root = "Roots", leaf = "Leaves", Flower = "Flowers")
tissue_annotation_color_map <- c(
  Root = "#CDB7FF",
  Leaf = "#FFE58F",
  Flower = "#FFC58F"
)

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
compartment_levels <- c("Cy", "Mt", "Pt")
compartment_color_map <- c(Cy = "#D55E00", Mt = "#0072B2", Pt = "#009E73")

# ComplexHeatmap 在"无设备"状态下测量文本宽度时，R 会临时打开默认 pdf 设备，而该设备的
# PostScript 字体库里没有 Arial，于是每个被测量的文本元素抛一条警告。测量值不受影响
# （实测与 cairo 下一致），且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。
# 只屏蔽这一条已知无害的警告，其他警告照常抛出。
muffle_arial_postscript_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("PostScript", msg) && grepl("Arial", msg)) invokeRestart("muffleWarning")
  })
}

# 文字统一走 shared 的无 showtext 路线：gpar 里显式指定 Arial
text_gp <- function(size, family = "Arial") {
  grid::gpar(fontsize = size, fontface = "bold", fontfamily = family, col = "black")
}

setup_heatmap_style <- function() {
  register_font()
}

read_structure_table <- function(rep_structure_path) {
  rep_structure <- readr::read_csv(rep_structure_path, show_col_types = FALSE)

  expected_cols <- c("ID", "Rep_ID", region_levels)
  missing_cols <- setdiff(expected_cols, names(rep_structure))
  if (length(missing_cols) > 0) {
    stop(
      "Representative tRNA structure CSV is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  rep_structure |>
    select(all_of(expected_cols)) |>
    mutate(
      display_Rep_ID = str_replace(Rep_ID, "^Nu-", "Cy-"),
      compartment = case_when(
        str_starts(Rep_ID, "Nu-") ~ "Cy",
        str_starts(Rep_ID, "Mt-") ~ "Mt",
        str_starts(Rep_ID, "Pt-") ~ "Pt",
        TRUE ~ NA_character_
      ),
      compartment = factor(compartment, levels = compartment_levels)
    )
}

# ============================================================================
# 输入读取与折叠（native 位 -> canonical 列）
# ============================================================================
read_union_sites <- function(file_path, value_col = default_value_col) {
  site_df <- readr::read_tsv(file_path, show_col_types = FALSE)

  required_cols <- c("chr", "pos", "fix_pos", value_col)
  missing_cols <- setdiff(required_cols, names(site_df))
  if (length(missing_cols) > 0) {
    stop(
      basename(file_path),
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  site_df |>
    transmute(
      chr,
      ID = chr,
      pos = as.integer(pos),
      canonical_pos = as.integer(fix_pos),
      value = as.numeric(.data[[value_col]]),
      tissues = if ("Tissues" %in% names(site_df)) as.character(site_df$Tissues) else NA_character_,
      peak_tissue = if ("peak_tissue" %in% names(site_df)) as.character(site_df$peak_tissue) else NA_character_
    ) |>
    filter(
      !is.na(ID),
      !is.na(canonical_pos),
      !is.na(value),
      canonical_pos >= 0L,
      canonical_pos <= max_canonical_pos
    )
}

report_non_representative <- function(site_df, structure_info) {
  dropped <- site_df |>
    anti_join(structure_info |> select(ID), by = "ID") |>
    distinct(ID) |>
    arrange(ID)

  if (nrow(dropped) > 0) {
    message(
      "Dropped ", nrow(dropped),
      " non-representative tRNA copies not present in All-tRNA_rep-structure.csv: ",
      paste(dropped$ID, collapse = ", ")
    )
  }

  invisible(dropped)
}

collapse_union_sites <- function(site_df, structure_info) {
  report_non_representative(site_df, structure_info)

  representative_df <- site_df |>
    inner_join(
      structure_info |> select(ID, Rep_ID, display_Rep_ID, compartment),
      by = "ID"
    )

  if (nrow(representative_df) == 0) {
    stop("No representative tRNA IDs remained after matching union `chr` values to the `ID` column.")
  }

  collapsed_df <- representative_df |>
    mutate(matrix_col_id = matrix_col_ids[canonical_pos + 1L]) |>
    group_by(ID, chr, Rep_ID, display_Rep_ID, compartment, matrix_col_id) |>
    summarise(
      value = max(value, na.rm = TRUE),
      n_sites = dplyr::n(),
      native_positions = paste(sort(unique(pos)), collapse = ","),
      tissues = paste(sort(unique(na.omit(tissues))), collapse = ";"),
      peak_tissues = paste(sort(unique(na.omit(peak_tissue))), collapse = ";"),
      canonical_pos = first(match(matrix_col_id, matrix_col_ids) - 1L),
      .groups = "drop"
    )

  collapsed_cells <- collapsed_df |> filter(n_sites > 1)
  if (nrow(collapsed_cells) > 0) {
    message(
      "Collapsed ", nrow(collapsed_cells),
      " heatmap cells with multiple union rows using max Delta Ratio."
    )
  }

  collapsed_df
}

build_row_info <- function(structure_info, collapsed_df) {
  structure_info |>
    filter(ID %in% collapsed_df$ID) |>
    mutate(compartment = factor(compartment, levels = compartment_levels)) |>
    arrange(compartment, match(ID, structure_info$ID))
}

build_union_matrix <- function(collapsed_df, row_info) {
  final_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_info),
    ncol = length(matrix_col_ids),
    dimnames = list(row_info$display_Rep_ID, matrix_col_ids)
  )

  if (nrow(collapsed_df) > 0) {
    row_idx <- match(collapsed_df$ID, row_info$ID)
    col_idx <- match(collapsed_df$matrix_col_id, matrix_col_ids)
    final_matrix[cbind(row_idx, col_idx)] <- collapsed_df$value
  }

  final_matrix
}

read_all_tissue_sites <- function(file_path, value_col = default_value_col) {
  site_df <- readr::read_tsv(file_path, show_col_types = FALSE)

  required_cols <- c("chr", "pos", "fix_pos", "Tissue", value_col)
  missing_cols <- setdiff(required_cols, names(site_df))
  if (length(missing_cols) > 0) {
    stop(
      basename(file_path),
      " is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  site_df |>
    transmute(
      Tissue = factor(Tissue, levels = tissue_levels, labels = tissue_display_map[tissue_levels]),
      chr,
      ID = chr,
      pos = as.integer(pos),
      canonical_pos = as.integer(fix_pos),
      value = as.numeric(.data[[value_col]])
    ) |>
    filter(
      !is.na(Tissue),
      !is.na(ID),
      !is.na(canonical_pos),
      !is.na(value),
      canonical_pos >= 0L,
      canonical_pos <= max_canonical_pos
    )
}

collapse_all_tissue_sites <- function(site_df, structure_info) {
  report_non_representative(site_df, structure_info)

  site_df |>
    inner_join(
      structure_info |> select(ID, Rep_ID, display_Rep_ID, compartment),
      by = "ID"
    ) |>
    mutate(matrix_col_id = matrix_col_ids[canonical_pos + 1L]) |>
    group_by(Tissue, ID, chr, Rep_ID, display_Rep_ID, compartment, matrix_col_id) |>
    summarise(
      value = max(value, na.rm = TRUE),
      n_sites = dplyr::n(),
      native_positions = paste(sort(unique(pos)), collapse = ","),
      canonical_pos = first(match(matrix_col_id, matrix_col_ids) - 1L),
      .groups = "drop"
    )
}

build_tissue_panel_row_info <- function(structure_info, collapsed_df) {
  collapsed_df |>
    distinct(Tissue, ID, Rep_ID, display_Rep_ID, compartment) |>
    mutate(
      Tissue = factor(Tissue, levels = unname(tissue_display_map[tissue_levels])),
      compartment = factor(compartment, levels = compartment_levels)
    ) |>
    arrange(Tissue, compartment, match(ID, structure_info$ID))
}

build_tissue_panel_matrix <- function(collapsed_df, row_info) {
  final_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_info),
    ncol = length(matrix_col_ids),
    dimnames = list(row_info$display_Rep_ID, matrix_col_ids)
  )

  if (nrow(collapsed_df) > 0) {
    row_keys <- paste(row_info$Tissue, row_info$ID, sep = "||")
    value_keys <- paste(collapsed_df$Tissue, collapsed_df$ID, sep = "||")
    row_idx <- match(value_keys, row_keys)
    col_idx <- match(collapsed_df$matrix_col_id, matrix_col_ids)
    final_matrix[cbind(row_idx, col_idx)] <- collapsed_df$value
  }

  final_matrix
}

# ============================================================================
# 注释（底部区段带 / 左侧分室 / 右侧组织）
# ============================================================================
build_bottom_annotation <- function() {
  tokens <- SIZE_TOKENS
  ComplexHeatmap::HeatmapAnnotation(
    Structure_Colors = ComplexHeatmap::anno_simple(
      x = as.character(column_metadata$region),
      col = region_color_map,
      height = grid::unit(3, "mm"),
      border = TRUE,
      gp = grid::gpar(col = "black", lwd = tokens$base_linewidth)
    ),
    Numbers = ComplexHeatmap::anno_text(
      x = column_metadata$display_number,
      rot = 90,
      just = "right",
      gp = text_gp(tokens$annotation_text_size)
    ),
    Regions = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = "white", col = NA),
      labels = layout_config$broad_label,
      labels_gp = text_gp(tokens$annotation_text_size),
      labels_rot = 0
    ),
    gap = grid::unit(1, "mm"),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
}

build_row_annotation <- function(row_info) {
  tokens <- SIZE_TOKENS
  ComplexHeatmap::rowAnnotation(
    Compartment = ComplexHeatmap::anno_simple(
      x = as.character(row_info$compartment),
      col = compartment_color_map,
      width = grid::unit(3.5, "mm"),
      border = TRUE,
      gp = grid::gpar(col = "black", lwd = tokens$base_linewidth)
    ),
    show_annotation_name = FALSE
  )
}

build_tissue_right_annotation <- function(tissue_split) {
  tokens <- SIZE_TOKENS
  tissue_labels <- levels(tissue_split)

  ComplexHeatmap::rowAnnotation(
    Tissue = ComplexHeatmap::anno_block(
      which = "row",
      labels = tissue_labels,
      gp = grid::gpar(
        fill = tissue_annotation_color_map[tissue_labels],
        col = "black",
        lwd = tokens$base_linewidth
      ),
      labels_gp = text_gp(tokens$annotation_text_size),
      labels_rot = 0,
      width = grid::unit(16, "mm")
    ),
    show_annotation_name = FALSE
  )
}

value_color_fun <- function() {
  circlize::colorRamp2(c(0, 1), c("#FFF5EE", "#8B0000"))
}

# ============================================================================
# 热图构造
# ============================================================================
build_union_heatmap <- function(final_matrix, row_info) {
  split_vec <- factor(column_metadata$region, levels = layout_config$region)
  compartment_split <- droplevels(row_info$compartment)

  old_dimname_padding <- ComplexHeatmap::ht_opt("DIMNAME_PADDING")
  on.exit(ComplexHeatmap::ht_opt(DIMNAME_PADDING = old_dimname_padding), add = TRUE)
  ComplexHeatmap::ht_opt(DIMNAME_PADDING = dimname_padding)

  ComplexHeatmap::Heatmap(
    final_matrix,
    name = default_value_title,
    col = value_color_fun(),
    na_col = "#FCFCFC",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = compartment_split,
    column_split = split_vec,
    column_gap = grid::unit(0, "mm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = SIZE_TOKENS$base_linewidth),
    left_annotation = build_row_annotation(row_info),
    bottom_annotation = build_bottom_annotation(),
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_names_gp = text_gp(SIZE_TOKENS$annotation_text_size),
    row_names_max_width = grid::unit(7, "cm"),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      title_gp = text_gp(SIZE_TOKENS$legend_title_size),
      labels_gp = text_gp(SIZE_TOKENS$legend_text_size)
    )
  )
}

build_tissue_panel_heatmap <- function(final_matrix, row_info) {
  split_vec <- factor(column_metadata$region, levels = layout_config$region)
  tissue_split <- droplevels(row_info$Tissue)

  old_dimname_padding <- ComplexHeatmap::ht_opt("DIMNAME_PADDING")
  on.exit(ComplexHeatmap::ht_opt(DIMNAME_PADDING = old_dimname_padding), add = TRUE)
  ComplexHeatmap::ht_opt(DIMNAME_PADDING = dimname_padding)

  ComplexHeatmap::Heatmap(
    final_matrix,
    name = default_value_title,
    col = value_color_fun(),
    na_col = "#FCFCFC",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = tissue_split,
    row_gap = grid::unit(2, "mm"),
    column_split = split_vec,
    column_gap = grid::unit(0, "mm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = SIZE_TOKENS$base_linewidth),
    left_annotation = build_row_annotation(row_info),
    right_annotation = build_tissue_right_annotation(tissue_split),
    bottom_annotation = build_bottom_annotation(),
    show_column_names = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_names_gp = text_gp(SIZE_TOKENS$annotation_text_size),
    row_names_max_width = grid::unit(7, "cm"),
    row_title = NULL,
    column_title = NULL,
    heatmap_legend_param = list(
      at = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
      title_gp = text_gp(SIZE_TOKENS$legend_title_size),
      labels_gp = text_gp(SIZE_TOKENS$legend_text_size)
    )
  )
}

# ============================================================================
# 绘制：热图铺满主视口，右侧独立视口放带边框的竖直色标
# ============================================================================
build_value_legend <- function(title = default_value_title) {
  ComplexHeatmap::Legend(
    col_fun = value_color_fun(),
    title = title,
    at = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0.00", "0.25", "0.50", "0.75", "1.00"),
    direction = "vertical",
    border = "black",
    legend_gp = grid::gpar(col = "black", lwd = SIZE_TOKENS$base_linewidth),
    background = "white",
    title_gp = text_gp(SIZE_TOKENS$legend_title_size),
    labels_gp = text_gp(SIZE_TOKENS$legend_text_size),
    grid_width = grid::unit(6, "mm"),
    legend_height = grid::unit(42, "mm"),
    tick_length = grid::unit(3, "mm")
  )
}

draw_heatmap_with_legend <- function(ht, legend_title = default_value_title) {
  legend_obj <- build_value_legend(legend_title)

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

save_heatmap <- function(ht, output_prefix, out_dir, width_in, height_in,
                         legend_title = default_value_title) {
  stem <- file.path(out_dir, paste0(output_prefix, "_heatmap"))

  muffle_arial_postscript_warnings(
    save_plot(
      plot = function() {
        grid::grid.newpage()
        draw_heatmap_with_legend(ht, legend_title)
      },
      filename = stem,
      width = width_in,
      height = height_in,
      formats = c("pdf", "png", "svg"),
      dpi = export_dpi
    )
  )

  list(
    plot_pdf = paste0(stem, ".pdf"),
    plot_png = paste0(stem, ".png"),
    plot_svg = paste0(stem, ".svg"),
    width_in = width_in,
    height_in = height_in
  )
}

# ============================================================================
# 对外的两个入口
# ============================================================================
require_input <- function(input_path, data_script) {
  if (!file.exists(input_path)) {
    stop(
      "Missing input table: ", input_path,
      "\n先跑 data_prep.sh 生成。",
      if (!is.null(data_script)) paste0("（", basename(data_script), "）") else ""
    )
  }
}

run_tissue_split_heatmap <- function(
  genotype,
  base_dir,
  out_dir,
  input_path = NULL,
  save_outputs = TRUE
) {
  setup_heatmap_style()

  rep_structure_path <- file.path(base_dir, "All-tRNA_rep-structure.csv")
  if (!file.exists(rep_structure_path)) {
    stop("Missing representative structure CSV: ", rep_structure_path)
  }

  if (is.null(input_path)) {
    input_path <- file.path(out_dir, paste0(genotype, "_all_tissues.dependent.tsv"))
  }
  require_input(input_path, file.path(out_dir, "data_prep.sh"))

  structure_info <- read_structure_table(rep_structure_path)
  all_tissue_sites <- read_all_tissue_sites(input_path)
  if (nrow(all_tissue_sites) == 0) {
    stop("No tissue-resolved dependent sites found in ", input_path)
  }

  collapsed_sites <- collapse_all_tissue_sites(all_tissue_sites, structure_info)
  if (nrow(collapsed_sites) == 0) {
    stop("No representative Cy/Mt/Pt rows remained after tissue-specific mapping.")
  }

  row_info <- build_tissue_panel_row_info(structure_info, collapsed_sites)
  final_matrix <- build_tissue_panel_matrix(collapsed_sites, row_info)
  heatmap_obj <- muffle_arial_postscript_warnings(
    build_tissue_panel_heatmap(final_matrix, row_info)
  )

  output_prefix <- paste(genotype, "delta_ratio", "tissue_split_dependent_trna", sep = "_")
  mapped_tsv <- file.path(out_dir, paste0(output_prefix, "_mapped_sites.tsv"))
  output_paths <- list(mapped_tsv = mapped_tsv)

  if (save_outputs) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(collapsed_sites, mapped_tsv)
    output_paths <- c(
      output_paths,
      save_heatmap(
        ht = heatmap_obj,
        output_prefix = output_prefix,
        out_dir = out_dir,
        width_in = 18,
        height_in = max(8.8, 0.19 * nrow(row_info) + 2.8)
      )
    )
  }

  metadata <- tibble::tibble(
    genotype = genotype,
    input_path = input_path,
    mapped_tsv = output_paths$mapped_tsv,
    plot_pdf = if (!is.null(output_paths$plot_pdf)) output_paths$plot_pdf else NA_character_,
    plot_png = if (!is.null(output_paths$plot_png)) output_paths$plot_png else NA_character_,
    plot_svg = if (!is.null(output_paths$plot_svg)) output_paths$plot_svg else NA_character_,
    n_input_rows = nrow(all_tissue_sites),
    n_rows = nrow(row_info)
  )

  list(
    genotype = genotype,
    out_dir = out_dir,
    metadata = metadata,
    input_sites = all_tissue_sites,
    mapped_sites = collapsed_sites,
    heatmap = heatmap_obj
  )
}

run_union_dependent_heatmap <- function(
  genotype,
  base_dir,
  out_dir,
  input_path = NULL,
  save_outputs = TRUE
) {
  setup_heatmap_style()

  rep_structure_path <- file.path(base_dir, "All-tRNA_rep-structure.csv")
  if (!file.exists(rep_structure_path)) {
    stop("Missing representative structure CSV: ", rep_structure_path)
  }

  if (is.null(input_path)) {
    input_path <- file.path(out_dir, paste0(genotype, "_union.dependent.tsv"))
  }
  require_input(input_path, file.path(out_dir, "data_prep.sh"))

  structure_info <- read_structure_table(rep_structure_path)
  union_sites <- read_union_sites(input_path)
  if (nrow(union_sites) == 0) {
    stop("No union dependent sites found in ", input_path)
  }

  collapsed_sites <- collapse_union_sites(union_sites, structure_info)

  if (nrow(collapsed_sites) == 0) {
    stop("No representative Cy/Mt/Pt tRNAs remained after filtering.")
  }

  row_info <- build_row_info(structure_info, collapsed_sites)
  final_matrix <- build_union_matrix(collapsed_sites, row_info)
  heatmap_obj <- muffle_arial_postscript_warnings(
    build_union_heatmap(final_matrix, row_info)
  )

  output_prefix <- paste(genotype, "delta_ratio", "union_dependent_trna", sep = "_")
  mapped_tsv <- file.path(out_dir, paste0(output_prefix, "_mapped_sites.tsv"))
  output_paths <- list(mapped_tsv = mapped_tsv)

  if (save_outputs) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(collapsed_sites, mapped_tsv)
    output_paths <- c(
      output_paths,
      save_heatmap(
        ht = heatmap_obj,
        output_prefix = output_prefix,
        out_dir = out_dir,
        width_in = 18,
        height_in = max(5.5, 0.19 * nrow(row_info) + 2.8)
      )
    )
  }

  metadata <- tibble::tibble(
    genotype = genotype,
    input_path = input_path,
    mapped_tsv = output_paths$mapped_tsv,
    plot_pdf = if (!is.null(output_paths$plot_pdf)) output_paths$plot_pdf else NA_character_,
    plot_png = if (!is.null(output_paths$plot_png)) output_paths$plot_png else NA_character_,
    plot_svg = if (!is.null(output_paths$plot_svg)) output_paths$plot_svg else NA_character_,
    n_input_rows = nrow(union_sites),
    n_tRNAs = nrow(row_info)
  )

  list(
    genotype = genotype,
    out_dir = out_dir,
    metadata = metadata,
    input_sites = union_sites,
    mapped_sites = collapsed_sites,
    heatmap = heatmap_obj,
    column_metadata = column_metadata,
    layout_config = layout_config
  )
}
