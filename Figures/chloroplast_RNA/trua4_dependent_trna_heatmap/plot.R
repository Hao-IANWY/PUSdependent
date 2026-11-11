#!/usr/bin/env Rscript
# trua4_dependent_trna_heatmap — structure heatmap of trua4-dependent
# chloroplast tRNA sites in leaf
#
# Inputs (beta-binomial calling criteria):
#   tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4-leaf.dependent.Sprinzl.tsv
#   Per-allele Mut_Ratio from trua4_{1,2}-leaf.dependent-sites.tsv (kept as
#   record columns only)
#
# Dependent sites are taken directly from the upstream official dependent
# table (AND across both mutant alleles, matching
# Figures/dependent_betabinomial); only ^tRNA-Pt- entries are kept. Canonical
# (Sprinzl) positions come from the Sprinzl_pos column of that table.
#
# Column layout = canonical positions 0-76 (17 regions from p1 through CCA,
# 77 columns), mirroring the heatmap in Figures/dependent_betabinomial/tRNA;
# rows = one dependent plastid tRNA per row.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

export_dpi <- 600
heatmap_draw_padding <- grid::unit(c(4, 16, 4, 54), "mm")

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

# When text width is measured with no open device, ComplexHeatmap temporarily
# opens the default pdf device, whose PostScript font library lacks Arial and
# emits one warning per measured text element. The measurements are unaffected;
# muffle only this known-benign warning.
muffle_arial_postscript_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    msg <- conditionMessage(w)
    if (grepl("PostScript", msg) && grepl("Arial", msg)) invokeRestart("muffleWarning")
  })
}

text_gp <- function(size) {
  grid::gpar(fontsize = size, fontface = "bold", fontfamily = "Arial", col = "black")
}

read_dependent_sites <- function(file_path) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c(
    "chr", "pos", "Sprinzl_pos", "strand", "context",
    "WT_Ratio", "Mut_Ratio", "Delta_Ratio(WT-Mut)"
  )
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      Sprinzl_pos = as.integer(Sprinzl_pos),
      strand = as.character(strand),
      context = as.character(context),
      WT_Ratio = as.numeric(WT_Ratio),
      Mut_Ratio = as.numeric(Mut_Ratio),
      `Delta_Ratio(WT-Mut)` = as.numeric(.data[["Delta_Ratio(WT-Mut)"]])
    ) %>%
    filter(grepl("^tRNA-Pt-", chr)) %>%
    transmute(
      chr,
      pos,
      Sprinzl_pos,
      strand,
      context,
      WT_Ratio,
      Mut_Ratio,
      `Delta_Ratio(WT-Mut)`,
      site_id = paste(chr, pos, strand, sep = ":")
    )
}

read_line_sites <- function(file_path, line_name) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "Mut_Ratio", "is_dependent")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand)
    ) %>%
    filter(grepl("^tRNA-Pt-", chr)) %>%
    transmute(
      site_id = paste(chr, pos, strand, sep = ":"),
      !!line_name := as.numeric(Mut_Ratio),
      !!paste0(line_name, "_dep") := as.logical(is_dependent)
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

build_value_legend <- function(col_fun, breaks, labels) {
  legend_height_mm <- max(38, 8.5 * length(breaks) + 6)

  ComplexHeatmap::Legend(
    col_fun = col_fun,
    title = "Delta Ratio",
    at = breaks,
    labels = labels,
    direction = "vertical",
    border = "black",
    legend_gp = grid::gpar(col = "black", lwd = SIZE_TOKENS$base_linewidth),
    background = "white",
    title_gp = text_gp(SIZE_TOKENS$legend_title_size),
    labels_gp = text_gp(SIZE_TOKENS$legend_text_size),
    grid_width = grid::unit(6, "mm"),
    legend_height = grid::unit(legend_height_mm, "mm"),
    tick_length = grid::unit(3, "mm")
  )
}

draw_heatmap_with_boxed_legend <- function(ht, col_fun, breaks, labels) {
  legend_obj <- build_value_legend(col_fun, breaks, labels)
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

  dependent_dir <- file.path(
    repo_root, "tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4"
  )
  dependent_path <- file.path(dependent_dir, "trua4-leaf.dependent.Sprinzl.tsv")
  line1_path <- file.path(dependent_dir, "trua4_1-leaf.dependent-sites.tsv")
  line2_path <- file.path(dependent_dir, "trua4_2-leaf.dependent-sites.tsv")

  required_paths <- c(dependent_path, line1_path, line2_path)
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0) {
    stop("Missing required inputs:\n", paste(missing_paths, collapse = "\n"))
  }

  dependent_tbl <- read_dependent_sites(dependent_path)
  if (nrow(dependent_tbl) == 0) {
    stop("No chloroplast (^tRNA-Pt-) sites in ", dependent_path)
  }

  union_tbl <- dependent_tbl %>%
    left_join(read_line_sites(line1_path, "trua4_1"), by = "site_id") %>%
    left_join(read_line_sites(line2_path, "trua4_2"), by = "site_id") %>%
    mutate(
      trua4_1_dep = coalesce(trua4_1_dep, FALSE),
      trua4_2_dep = coalesce(trua4_2_dep, FALSE),
      dependency_support = case_when(
        trua4_1_dep & trua4_2_dep ~ "shared",
        trua4_1_dep ~ "trua4_1_only",
        trua4_2_dep ~ "trua4_2_only",
        TRUE ~ "none"
      )
    ) %>%
    arrange(Sprinzl_pos, chr, pos) %>%
    select(
      chr, pos, Sprinzl_pos, strand, context,
      WT_Ratio, Mut_Ratio, `Delta_Ratio(WT-Mut)`,
      trua4_1, trua4_2, dependency_support
    )

  readr::write_tsv(
    union_tbl,
    file.path(script_dir, "trua4_pt_trna_dependent_union.tsv")
  )

  mapped_tbl <- union_tbl %>%
    mutate(
      display_label = str_replace(chr, "^tRNA-", ""),
      canonical_pos = as.integer(Sprinzl_pos),
      matrix_col_id = matrix_col_ids[canonical_pos + 1L]
    ) %>%
    filter(!is.na(matrix_col_id), canonical_pos >= 0L, canonical_pos <= max_canonical_pos) %>%
    arrange(canonical_pos, display_label)

  if (nrow(mapped_tbl) == 0) {
    stop("No chloroplast tRNA rows mapped onto the canonical 0-76 layout.")
  }

  readr::write_tsv(
    mapped_tbl,
    file.path(script_dir, "trua4_pt_trna_dependent_mapped_sites.tsv")
  )

  row_info <- mapped_tbl %>%
    distinct(chr, display_label, canonical_pos) %>%
    arrange(canonical_pos, display_label)

  heatmap_matrix <- matrix(
    NA_real_,
    nrow = nrow(row_info),
    ncol = length(matrix_col_ids),
    dimnames = list(row_info$display_label, matrix_col_ids)
  )

  for (i in seq_len(nrow(mapped_tbl))) {
    row_idx <- match(mapped_tbl$display_label[[i]], row_info$display_label)
    col_idx <- match(mapped_tbl$matrix_col_id[[i]], matrix_col_ids)
    current_value <- heatmap_matrix[row_idx, col_idx]
    new_value <- mapped_tbl[["Delta_Ratio(WT-Mut)"]][[i]]
    if (is.na(current_value) || new_value > current_value) {
      heatmap_matrix[row_idx, col_idx] <- new_value
    }
  }

  value_axis <- y_axis_cover(max(mapped_tbl[["Delta_Ratio(WT-Mut)"]], na.rm = TRUE), n = 4)
  col_fun <- circlize::colorRamp2(c(0, value_axis$top), c("#FFF5EE", "#8B0000"))
  legend_labels <- formatC(value_axis$breaks, format = "f", digits = 2)

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
        at = value_axis$breaks,
        labels = legend_labels,
        title_gp = text_gp(SIZE_TOKENS$legend_title_size),
        labels_gp = text_gp(SIZE_TOKENS$legend_text_size)
      )
    )
  )

  width_in <- 18
  height_in <- max(6.0, 0.55 * nrow(row_info) + 2.8)
  stem <- file.path(script_dir, "trua4_pt_trna_dependent_heatmap")

  muffle_arial_postscript_warnings(
    save_plot(
      plot = function() {
        grid::grid.newpage()
        draw_heatmap_with_boxed_legend(
          ht,
          col_fun = col_fun,
          breaks = value_axis$breaks,
          labels = legend_labels
        )
      },
      filename = stem,
      width = width_in,
      height = height_in,
      formats = c("pdf", "png", "svg"),
      dpi = export_dpi
    )
  )
}

if (sys.nframe() == 0) {
  run_analysis()
}
