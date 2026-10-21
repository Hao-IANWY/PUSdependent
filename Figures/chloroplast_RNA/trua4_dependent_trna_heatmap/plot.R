#!/usr/bin/env Rscript
# trua4_dependent_trna_heatmap — trua4 dependent 叶绿体 tRNA Ψ 位点的结构热图
#
# 移植自 archive/20260902/paper/chloroplast_RNA/trua4_dependent_trna_heatmap/plot.R。
# 输入换成当前 methylKit 口径：
#   tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv
#   （只取 is_dependent == TRUE 且 chr 以 tRNA-Pt- 开头的位点）
# 结构布局沿用 All-tRNA_rep-structure.csv（46 Nu / 17 Mt / 30 Pt）。
#
# 与 archive 版的两处差异：
#   1. 当前 dependent-sites 表没有 fix_pos，按结构表区段长度现算 native pos -> canonical 0-76。
#   2. 色标上限按数据自适应（至少 0.3），避免高 Delta 位点全部饱和在同一个颜色。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  library(tibble)
  library(grid)
})

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

column_metadata <- layout_config %>%
  mutate(rel_pos = lapply(fixed_length, seq_len)) %>%
  tidyr::unnest(rel_pos) %>%
  mutate(
    matrix_col_id = paste(region, rel_pos, sep = "_"),
    display_number = 0:(dplyr::n() - 1L)
  )

matrix_col_ids <- column_metadata$matrix_col_id
max_canonical_pos <- max(column_metadata$display_number)
region_offset <- setNames(cumsum(c(0L, head(layout_config$fixed_length, -1L))), region_levels)
region_color_map <- setNames(layout_config$color_hex, as.character(layout_config$region))

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]
  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }
  normalizePath(getwd())
}

# ComplexHeatmap 在"无设备"状态测量文本宽度时会临时打开默认 pdf 设备，
# 该设备字体库没有 Arial；测量值不受影响，只屏蔽这一条已知无害警告。
muffle_arial_postscript_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("PostScript", msg) && grepl("Arial", msg)) invokeRestart("muffleWarning")
  })
}

text_gp <- function(size) {
  grid::gpar(fontsize = size, fontface = "bold", fontfamily = "Arial", col = "black")
}

# native pos -> canonical 0-76：按结构表区段顺序，把每个区段的实际长度
# 重标定到固定长度（与 build_genotype_dependent_tables.sh / archive 的 awk 一致，
# 半值向上取整）。
build_canonical_map <- function(rep_structure) {
  rows <- vector("list", nrow(rep_structure))
  for (i in seq_len(nrow(rep_structure))) {
    rep_id <- rep_structure$Rep_ID[[i]]
    native_start <- 1L
    pieces <- list()
    for (region in region_levels) {
      seq <- rep_structure[[region]][[i]]
      region_len <- nchar(seq)
      if (is.na(region_len) || region_len == 0) next
      fixed_len <- layout_config$fixed_length[[match(region, region_levels)]]
      native_pos <- native_start:(native_start + region_len - 1L)
      if (fixed_len == 1L || region_len == 1L) {
        canonical_region_pos <- rep(1L, region_len)
      } else {
        canonical_region_pos <- as.integer(
          floor((seq_len(region_len) - 1L) * (fixed_len - 1L) / (region_len - 1L) + 0.5)
        ) + 1L
      }
      pieces[[length(pieces) + 1L]] <- tibble(
        Rep_ID = rep_id,
        native_pos = native_pos,
        canonical_pos = region_offset[[region]] + canonical_region_pos - 1L
      )
      native_start <- native_start + region_len
    }
    rows[[i]] <- bind_rows(pieces)
  }
  bind_rows(rows) %>% distinct(Rep_ID, native_pos, .keep_all = TRUE)
}

read_dependent_line <- function(file_path, line_name) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  delta_col <- "Delta_Ratio(WT-Mut)"
  required_cols <- c("chr", "pos", "strand", "context", "WT_Ratio", "Mut_Ratio", delta_col, "is_dependent")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      context = as.character(context),
      WT_Ratio = as.numeric(WT_Ratio),
      Mut_Ratio = as.numeric(Mut_Ratio),
      delta_ratio = as.numeric(.data[[delta_col]]),
      is_dependent = as.logical(is_dependent)
    ) %>%
    filter(grepl("^tRNA-Pt-", chr), !is.na(is_dependent), is_dependent) %>%
    transmute(
      chr,
      pos,
      strand,
      context,
      WT_Ratio,
      !!line_name := Mut_Ratio,
      !!paste0(line_name, "_delta") := delta_ratio
    )
}

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
      # 77 个位号挤在 18 in 宽里；16pt 会连成一条黑带，12pt 才逐位可读。
      gp = text_gp(SIZE_TOKENS$annotation_text_size - 4)
    ),
    Regions = ComplexHeatmap::anno_block(
      gp = grid::gpar(fill = "white", col = NA),
      labels = layout_config$broad_label,
      # 区段名比位号小 2pt：Anticodon Stem/Loop/Stem 三个窄块并排时
      # 16pt 的 "Anticodon" 会互相压字，14pt 才能在每个块宽内放下。
      labels_gp = text_gp(SIZE_TOKENS$annotation_text_size - 2),
      labels_rot = 0
    ),
    gap = grid::unit(1, "mm"),
    show_legend = FALSE,
    show_annotation_name = FALSE
  )
}

build_value_legend <- function(ramp_max, legend_breaks) {
  legend_obj <- ComplexHeatmap::Legend(
    col_fun = circlize::colorRamp2(c(0, ramp_max), c("#FFF5EE", "#8B0000")),
    title = "Delta Ratio",
    at = legend_breaks,
    labels = formatC(legend_breaks, format = "f", digits = 2),
    direction = "vertical",
    border = "black",
    legend_gp = grid::gpar(col = "black", lwd = SIZE_TOKENS$base_linewidth),
    background = "white",
    title_gp = text_gp(SIZE_TOKENS$legend_title_size),
    labels_gp = text_gp(SIZE_TOKENS$legend_text_size),
    grid_width = grid::unit(6, "mm"),
    legend_height = grid::unit(max(42, 9 * length(legend_breaks)), "mm"),
    tick_length = grid::unit(3, "mm")
  )
  legend_obj
}

draw_heatmap_with_boxed_legend <- function(ht, ramp_max, legend_breaks) {
  legend_obj <- build_value_legend(ramp_max, legend_breaks)
  heatmap_draw_padding <- grid::unit(c(4, 16, 4, 54), "mm")
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

run_analysis <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  register_font()

  line1_path <- file.path(repo_root, "tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_1-leaf.dependent-sites.tsv")
  line2_path <- file.path(repo_root, "tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_2-leaf.dependent-sites.tsv")
  structure_path <- file.path(script_dir, "All-tRNA_rep-structure.csv")

  required_paths <- c(line1_path, line2_path, structure_path)
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0) {
    stop("Missing required inputs:\n", paste(missing_paths, collapse = "\n"))
  }

  structure_info <- readr::read_csv(structure_path, show_col_types = FALSE) %>%
    select(ID, Rep_ID, all_of(region_levels)) %>%
    filter(stringr::str_starts(Rep_ID, "Pt-")) %>%
    mutate(display_label = stringr::str_replace(ID, "^tRNA-", ""))

  canonical_map <- build_canonical_map(structure_info)

  line1_tbl <- read_dependent_line(line1_path, "trua4_1")
  line2_tbl <- read_dependent_line(line2_path, "trua4_2")

  union_tbl <- full_join(
    line1_tbl,
    line2_tbl,
    by = c("chr", "pos", "strand"),
    suffix = c("_x", "_y")
  ) %>%
    mutate(
      context = coalesce(context_x, context_y),
      WT_Ratio = coalesce(WT_Ratio_x, WT_Ratio_y),
      trua4_1 = as.numeric(trua4_1),
      trua4_2 = as.numeric(trua4_2),
      Mut_Ratio = case_when(
        !is.na(trua4_1) & !is.na(trua4_2) ~ (trua4_1 + trua4_2) / 2,
        !is.na(trua4_1) ~ trua4_1,
        !is.na(trua4_2) ~ trua4_2,
        TRUE ~ NA_real_
      ),
      dependency_support = case_when(
        !is.na(trua4_1) & !is.na(trua4_2) ~ "shared",
        !is.na(trua4_1) ~ "trua4_1_only",
        !is.na(trua4_2) ~ "trua4_2_only",
        TRUE ~ "none"
      ),
      `Delta_Ratio(WT-Mut)` = WT_Ratio - Mut_Ratio
    ) %>%
    filter(!is.na(Mut_Ratio)) %>%
    mutate(rep_id = str_remove(chr, "^tRNA-") %>% str_remove("-[0-9]+$")) %>%
    left_join(canonical_map, by = c("rep_id" = "Rep_ID", "pos" = "native_pos")) %>%
    group_by(chr, pos, strand) %>%
    slice_min(canonical_pos, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      chr,
      pos,
      fix_pos = canonical_pos,
      strand,
      context,
      WT_Ratio,
      Mut_Ratio,
      `Delta_Ratio(WT-Mut)`,
      trua4_1,
      trua4_2,
      dependency_support
    ) %>%
    arrange(chr, fix_pos, pos)

  if (nrow(union_tbl) == 0) {
    stop("No chloroplast trua4-dependent tRNA sites remained after union.")
  }
  if (anyNA(union_tbl$fix_pos)) {
    bad_rows <- filter(union_tbl, is.na(fix_pos))
    stop(
      "Failed to map native positions to the representative structure layout:\n",
      paste(unique(paste0(bad_rows$chr, ":", bad_rows$pos)), collapse = "\n")
    )
  }

  readr::write_tsv(union_tbl, file.path(script_dir, "trua4_pt_trna_dependent_union.tsv"))

  mapped_tbl <- union_tbl %>%
    inner_join(
      structure_info %>% select(ID, Rep_ID, display_label),
      by = c("chr" = "ID")
    ) %>%
    mutate(
      canonical_pos = fix_pos,
      matrix_col_id = matrix_col_ids[canonical_pos + 1L],
      dependency_support = factor(dependency_support, levels = c("shared", "trua4_1_only", "trua4_2_only"))
    ) %>%
    filter(!is.na(matrix_col_id), canonical_pos >= 0L, canonical_pos <= max_canonical_pos)

  if (nrow(mapped_tbl) == 0) {
    stop("No chloroplast tRNA rows matched representative Pt-tRNA structures.")
  }

  readr::write_tsv(mapped_tbl, file.path(script_dir, "trua4_pt_trna_dependent_mapped_sites.tsv"))

  row_info <- mapped_tbl %>%
    distinct(chr, display_label) %>%
    mutate(order_id = match(chr, structure_info$ID)) %>%
    arrange(order_id, chr)

  heatmap_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_info),
    ncol = length(matrix_col_ids),
    dimnames = list(row_info$display_label, matrix_col_ids)
  )

  for (i in seq_len(nrow(mapped_tbl))) {
    row_idx <- match(mapped_tbl$chr[[i]], row_info$chr)
    col_idx <- match(mapped_tbl$matrix_col_id[[i]], matrix_col_ids)
    current_value <- heatmap_matrix[row_idx, col_idx]
    new_value <- mapped_tbl$`Delta_Ratio(WT-Mut)`[[i]]
    if (is.na(current_value) || new_value > current_value) {
      heatmap_matrix[row_idx, col_idx] <- new_value
    }
  }

  ramp_max <- max(0.3, ceiling(max(mapped_tbl$`Delta_Ratio(WT-Mut)`, na.rm = TRUE) * 10) / 10)
  legend_breaks <- seq(0, ramp_max, by = 0.1)
  col_fun <- circlize::colorRamp2(c(0, ramp_max), c("#FFF5EE", "#8B0000"))

  ht <- muffle_arial_postscript_warnings(
    ComplexHeatmap::Heatmap(
      heatmap_matrix,
      name = "Delta Ratio",
      col = col_fun,
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
      row_names_max_width = grid::unit(9, "cm"),
      row_title = NULL,
      column_title = NULL,
      heatmap_legend_param = list(
        at = legend_breaks,
        labels = formatC(legend_breaks, format = "f", digits = 2)
      )
    )
  )

  height_in <- max(4.8, 0.55 * nrow(row_info) + 2.8)
  muffle_arial_postscript_warnings(
    save_plot(
      plot = function() {
        grid::grid.newpage()
        draw_heatmap_with_boxed_legend(ht, ramp_max, legend_breaks)
      },
      filename = file.path(script_dir, "trua4_pt_trna_dependent_heatmap"),
      width = 18,
      height = height_in,
      dpi = 600
    )
  )
}

if (sys.nframe() == 0) {
  run_analysis()
}
