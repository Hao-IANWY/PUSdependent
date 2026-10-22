#!/usr/bin/env Rscript
# plot.R — 跨基因型 merged dependent tRNA 修饰比例热图（Root / leaf / Flower 各一张）
#
# 移植自 archive/20260902/paper/tRNA_dependent/merge_heatmap/plot.R，样式改走 Figures/shared。
# 与 archive 的三处差异：
#   1. 位点来源换成本目录各基因型的 <genotype>_all_tissues.dependent.tsv（methylKit 口径，自带 fix_pos）
#   2. WT/突变体修饰比例全部回填自 tRNA/STAT/adjustGap/ncRNA 的 adjusted.tsv（data_prep.sh 聚合）；
#      archive 的 WT 用 Fisher.tsv，但 tRNA 流程的 Fisher.tsv 是过滤子集，不能全位点回填
#   3. 分室 Cy/Mt/Pt 全进图（2026-10-12 老板拍板；archive 只 Cy/Mt）
#
# 图形：行 = 代表 tRNA（按分室分组），列 = fix_pos × 样本（WT + 7 突变体），
#       格值 = 修饰比例（0-1，缺失浅灰），dependent 格标 "*"。
#
# 输入：
#   ../<genotype>/<genotype>_all_tissues.dependent.tsv   7 个基因型（位点并集 + dependent "*" 标记）
#   merged_sample_ratios.tsv                              data_prep.sh 产物
# 输出（每组织）：
#   <tissue>_merged_dependent_site_ratio_heatmap.{pdf,png,svg}
#   <tissue>_merged_dependent_site_ratio_long.tsv / _matrix.tsv

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
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", cmd_args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", cmd_args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
tRNA_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
shared_dir <- normalizePath(file.path(tRNA_dir, "..", "..", "shared"), mustWork = TRUE)

source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))
source(file.path(tRNA_dir, "dependent_trna_heatmap_common.R"))

genotypes <- c("trua1a", "trua1b", "trua2", "trua3", "trua4", "trua5", "trua6")
tissues <- c("Root", "leaf", "Flower")
sample_levels <- c("WT", genotypes)
ratio_input <- file.path(script_dir, "merged_sample_ratios.tsv")

# ============================================================================
# 输入读取与并集
# ============================================================================
read_genotype_tissue_sources <- function(genotype) {
  input_path <- file.path(tRNA_dir, genotype, paste0(genotype, "_all_tissues.dependent.tsv"))
  if (!file.exists(input_path)) {
    stop("Missing tissue-resolved dependent-site input: ", input_path,
         "\n先跑 Figures/dependent_methylKit/tRNA/", genotype, "/data_prep.sh")
  }

  site_df <- readr::read_tsv(input_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "fix_pos", "strand", "context", "WT_Ratio", "Mut_Ratio", "Tissue")
  missing_cols <- setdiff(required_cols, names(site_df))
  if (length(missing_cols) > 0) {
    stop(basename(input_path), " is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  site_df |>
    transmute(
      source_genotype = genotype,
      Tissue = as.character(Tissue),
      chr,
      ID = chr,
      pos = as.integer(pos),
      fix_pos = as.integer(fix_pos),
      strand = as.character(strand),
      context = as.character(context),
      WT_Ratio = as.numeric(WT_Ratio),
      Mut_Ratio = as.numeric(Mut_Ratio)
    ) |>
    filter(!is.na(Tissue), !is.na(ID), !is.na(fix_pos), !is.na(WT_Ratio), !is.na(Mut_Ratio))
}

build_union_site_list <- function(tissue_sources) {
  tissue_sources |>
    group_by(Tissue, ID, chr, pos, fix_pos, strand) |>
    summarise(
      source_genotypes = paste(sort(unique(source_genotype)), collapse = ","),
      contexts = paste(sort(unique(na.omit(context))), collapse = ";"),
      .groups = "drop"
    )
}

attach_structure_info <- function(collapsed_df, structure_info) {
  non_representative <- collapsed_df |>
    anti_join(structure_info |> select(ID), by = "ID") |>
    distinct(ID) |>
    arrange(ID)

  if (nrow(non_representative) > 0) {
    message(
      "Dropped ", nrow(non_representative),
      " non-representative tRNA IDs not present in All-tRNA_rep-structure.csv: ",
      paste(non_representative$ID, collapse = ", ")
    )
  }

  collapsed_df |>
    inner_join(
      structure_info |> select(ID, Rep_ID, display_Rep_ID, compartment),
      by = "ID"
    ) |>
    mutate(compartment = factor(compartment, levels = compartment_levels))
}

read_sample_ratios <- function(path) {
  if (!file.exists(path)) {
    stop("Missing sample ratio table: ", path, "\n先跑 data_prep.sh 生成。")
  }

  ratio_df <- readr::read_tsv(path, show_col_types = FALSE)
  required_cols <- c("tissue", "sample", "chr", "pos", "strand", "ratio")
  missing_cols <- setdiff(required_cols, names(ratio_df))
  if (length(missing_cols) > 0) {
    stop(basename(path), " is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  ratio_df |>
    transmute(
      Tissue = as.character(tissue),
      sample = as.character(sample),
      chr,
      pos = as.integer(pos),
      strand = as.character(strand),
      ratio = as.numeric(ratio)
    ) |>
    filter(sample %in% sample_levels)
}

# ============================================================================
# 位点 × 样本 网格与矩阵
# ============================================================================
build_mapped_ratio_table <- function(union_sites, tissue_ratios, structure_info, sample_levels) {
  union_sites |>
    select(Tissue, ID, chr, pos, fix_pos, strand, source_genotypes, contexts) |>
    tidyr::crossing(sample = sample_levels) |>
    left_join(tissue_ratios, by = c("chr", "pos", "strand", "sample")) |>
    inner_join(
      structure_info |> select(ID, Rep_ID, display_Rep_ID, compartment),
      by = "ID"
    ) |>
    mutate(compartment = factor(compartment, levels = compartment_levels)) |>
    group_by(Tissue, ID, chr, fix_pos, sample, Rep_ID, display_Rep_ID, compartment) |>
    summarise(
      ratio = if (all(is.na(ratio))) NA_real_ else max(ratio, na.rm = TRUE),
      native_positions = paste(sort(unique(pos)), collapse = ","),
      source_genotypes = paste(sort(unique(source_genotypes)), collapse = ","),
      contexts = paste(sort(unique(na.omit(contexts))), collapse = ";"),
      .groups = "drop"
    )
}

build_merged_long <- function(mapped_df, structure_info, sample_levels, tissue) {
  row_info <- structure_info |>
    filter(ID %in% mapped_df$ID) |>
    mutate(compartment = factor(compartment, levels = compartment_levels)) |>
    arrange(compartment, match(ID, structure_info$ID)) |>
    select(ID, Rep_ID, display_Rep_ID, compartment)

  fix_pos_levels <- sort(unique(mapped_df$fix_pos))

  complete_grid <- tidyr::expand_grid(
    ID = row_info$ID,
    fix_pos = fix_pos_levels,
    sample = sample_levels
  )

  complete_grid |>
    left_join(row_info, by = "ID") |>
    left_join(mapped_df, by = c("ID", "fix_pos", "sample", "Rep_ID", "display_Rep_ID", "compartment")) |>
    mutate(
      Tissue = tissue,
      sample = factor(sample, levels = sample_levels),
      column_id = paste0("fix", fix_pos, "__", sample)
    ) |>
    arrange(compartment, match(ID, row_info$ID), fix_pos, sample)
}

build_ratio_matrix <- function(merged_long, row_info, column_info) {
  ratio_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_info),
    ncol = nrow(column_info),
    dimnames = list(row_info$display_Rep_ID, column_info$column_id)
  )

  value_rows <- merged_long |>
    filter(!is.na(ratio))

  if (nrow(value_rows) > 0) {
    row_idx <- match(value_rows$ID, row_info$ID)
    col_idx <- match(value_rows$column_id, column_info$column_id)
    ratio_matrix[cbind(row_idx, col_idx)] <- value_rows$ratio
  }

  ratio_matrix
}

build_dependent_marker_matrix <- function(merged_long, row_info, column_info) {
  marker_matrix <- matrix(
    FALSE,
    nrow = nrow(row_info),
    ncol = nrow(column_info),
    dimnames = list(row_info$display_Rep_ID, column_info$column_id)
  )

  marker_rows <- merged_long |>
    filter(
      !is.na(source_genotypes),
      as.character(sample) != "WT",
      stringr::str_detect(
        source_genotypes,
        paste0("(^|,)", stringr::str_replace_all(as.character(sample), "([\\W])", "\\\\\\1"), "(,|$)")
      )
    )

  if (nrow(marker_rows) > 0) {
    row_idx <- match(marker_rows$ID, row_info$ID)
    col_idx <- match(marker_rows$column_id, column_info$column_id)
    marker_matrix[cbind(row_idx, col_idx)] <- TRUE
  }

  marker_matrix
}

# ============================================================================
# 热图构造
# ============================================================================
build_merged_heatmap <- function(ratio_matrix, marker_matrix, row_info, column_info) {
  column_split <- factor(column_info$fix_pos, levels = unique(column_info$fix_pos))
  column_position_labels <- paste0("tRNA_", unique(column_info$fix_pos))
  row_labels <- stringr::str_remove(rownames(ratio_matrix), "^(Cy|Mt|Pt)-")

  old_dimname_padding <- ComplexHeatmap::ht_opt("DIMNAME_PADDING")
  on.exit(ComplexHeatmap::ht_opt(DIMNAME_PADDING = old_dimname_padding), add = TRUE)
  ComplexHeatmap::ht_opt(DIMNAME_PADDING = dimname_padding)

  ComplexHeatmap::Heatmap(
    ratio_matrix,
    name = "Ratio",
    col = value_color_fun(),
    na_col = "#FCFCFC",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = droplevels(row_info$compartment),
    column_split = column_split,
    column_gap = grid::unit(0.7, "mm"),
    border = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = SIZE_TOKENS$base_linewidth),
    left_annotation = build_row_annotation(row_info),
    show_column_names = TRUE,
    column_labels = as.character(column_info$sample),
    column_names_rot = 90,
    # 样本列密（fix_pos × 8），字号无法套用 tokens 档位，沿用 archive 的 7pt
    column_names_gp = text_gp(7),
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_names_side = "left",
    row_labels = row_labels,
    row_names_gp = text_gp(SIZE_TOKENS$annotation_text_size),
    row_names_max_width = grid::unit(7, "cm"),
    row_title = NULL,
    column_title = column_position_labels,
    column_title_gp = text_gp(SIZE_TOKENS$annotation_text_size),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (isTRUE(marker_matrix[i, j])) {
        grid::grid.text(
          "*",
          x = x,
          y = y - height * 0.10,
          gp = grid::gpar(fontsize = 9, fontface = "bold", fontfamily = "Arial", col = "black")
        )
      }
    },
    show_heatmap_legend = FALSE
  )
}

# ============================================================================
# 主流程
# ============================================================================
invisible(setup_heatmap_style())

structure_info <- read_structure_table(file.path(tRNA_dir, "All-tRNA_rep-structure.csv"))
tissue_sources <- purrr::map_dfr(genotypes, read_genotype_tissue_sources)
union_sites <- build_union_site_list(tissue_sources)
mapped_union_sites <- attach_structure_info(union_sites, structure_info)

if (nrow(mapped_union_sites) == 0) {
  stop("No representative tRNA rows remained after merging dependent-site inputs.")
}

all_ratios <- read_sample_ratios(ratio_input)

metadata <- purrr::map_dfr(tissues, function(tissue) {
  tissue_union_sites <- mapped_union_sites |>
    filter(Tissue == tissue)

  if (nrow(tissue_union_sites) == 0) {
    warning("No mapped dependent sites for tissue: ", tissue)
    return(tibble::tibble(
      tissue = tissue,
      n_union_tRNA_sites = 0L,
      n_displayed_tRNAs = 0L,
      n_corrected_positions = 0L,
      n_heatmap_columns = 0L
    ))
  }

  tissue_ratios <- all_ratios |>
    filter(Tissue == tissue) |>
    select(-Tissue)

  tissue_df <- build_mapped_ratio_table(tissue_union_sites, tissue_ratios, structure_info, sample_levels)
  merged_long <- build_merged_long(tissue_df, structure_info, sample_levels, tissue)

  row_info <- structure_info |>
    filter(ID %in% tissue_df$ID) |>
    mutate(compartment = factor(compartment, levels = compartment_levels)) |>
    arrange(compartment, match(ID, structure_info$ID)) |>
    select(ID, Rep_ID, display_Rep_ID, compartment)

  column_info <- merged_long |>
    distinct(fix_pos, sample, column_id) |>
    mutate(sample = factor(sample, levels = sample_levels)) |>
    arrange(fix_pos, sample)

  ratio_matrix <- build_ratio_matrix(merged_long, row_info, column_info)
  marker_matrix <- build_dependent_marker_matrix(merged_long, row_info, column_info)
  heatmap_obj <- muffle_arial_postscript_warnings(
    build_merged_heatmap(ratio_matrix, marker_matrix, row_info, column_info)
  )

  output_prefix <- paste0(tolower(tissue), "_merged_dependent_site_ratio")
  n_union_tRNA_sites <- nrow(dplyr::distinct(tissue_union_sites, ID, fix_pos))

  readr::write_tsv(
    merged_long |>
      select(
        Tissue, ID, Rep_ID, display_Rep_ID, compartment, fix_pos, sample, ratio,
        native_positions, source_genotypes, contexts
      ),
    file.path(script_dir, paste0(output_prefix, "_long.tsv"))
  )
  readr::write_tsv(
    ratio_matrix |>
      as.data.frame(check.names = FALSE) |>
      tibble::rownames_to_column("display_Rep_ID"),
    file.path(script_dir, paste0(output_prefix, "_matrix.tsv"))
  )
  output_paths <- save_heatmap(
    ht = heatmap_obj,
    output_prefix = output_prefix,
    out_dir = script_dir,
    width_in = max(18, 0.22 * ncol(ratio_matrix) + 4),
    height_in = max(6, 0.12 * nrow(ratio_matrix) + 2.8),
    legend_title = "Ratio"
  )

  tibble::tibble(
    tissue = tissue,
    n_union_tRNA_sites = n_union_tRNA_sites,
    n_displayed_tRNAs = nrow(row_info),
    n_corrected_positions = dplyr::n_distinct(column_info$fix_pos),
    n_heatmap_columns = ncol(ratio_matrix),
    plot_pdf = output_paths$plot_pdf,
    plot_png = output_paths$plot_png,
    plot_svg = output_paths$plot_svg
  )
})

print(metadata)
message("[plot] Done: merged dependent tRNA ratio heatmaps (Root / leaf / Flower)")
