#!/usr/bin/env Rscript
# site_specificity_common.R — 位点特异性（site specificity）图组公共实现
#
# 输入是 build_site_specificity_table.sh 生成的宽表（一行一个 WT 修饰的
# 代表 tRNA 位点，自带每条株系的修饰比与 per-line dependent 判定）：
#   chr, pos, fix_pos, strand, context, WT_Ratio, Treat_avg_depth,
#   对每个基因型 g：<g>_1, <g>_2, <g>_1_dependent, <g>_2_dependent
#
# 两张图：
#   1. 热图（样式对齐 ../mRNA/heatmap/plot.R）：行 = 入选 tRNA（WT 比例降序），
#      列 = WT + 各突变体株系，格值 = 修饰比例 0–100%，红框 = 该株系
#      beta-binomial per-line 判定为 dependent；缺失格子灰底 NA。
#   2. 修饰比例分布（样式对齐 ../dependent_trna_distribution_common.R 的
#      箱线 + 抖动点）：x = WT + 各株系，y = 修饰比例(%)。
#
# 只被各组织目录的 plot.R source，不直接运行。

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(ggtext)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# 每基因型两条株系的一对填充色（浅→深），WT 灰与基因型图组同款。
# 新基因型进图时在此补一对，缺色直接报错。
site_specificity_line_palette <- list(
  trua1a = c("#F4A3AECC", "#E26D7FCC"),
  trua5 = c("#A3C8E8CC", "#3E7CB1CC"),
  trua3 = c("#8FD3B0CC", "#009E73CC"),
  trua6 = c("#C9B3E8CC", "#7B52ABCC")
)

site_specificity_line_cols <- function(genotypes) {
  unlist(lapply(genotypes, function(g) paste0(g, "_", c("1", "2"))),
         use.names = FALSE)
}

site_specificity_fill_map <- function(genotypes) {
  fill <- c(WT = "#9E9E9ECC")
  for (g in genotypes) {
    pair <- site_specificity_line_palette[[g]]
    if (is.null(pair)) {
      stop("No line palette for genotype '", g,
           "'; add one pair to site_specificity_line_palette.")
    }
    fill[[paste0(g, "_1")]] <- pair[[1]]
    fill[[paste0(g, "_2")]] <- pair[[2]]
  }
  fill
}

prepare_site_specificity_table <- function(path, genotypes) {
  line_cols <- site_specificity_line_cols(genotypes)
  dep_cols <- paste0(line_cols, "_dependent")
  required_cols <- c("chr", "pos", "fix_pos", "strand", "context", "WT_Ratio",
                     "Treat_avg_depth", line_cols, dep_cols)
  if (!file.exists(path)) {
    stop("Input table not found: ", path, "\n先跑 data_prep.sh 生成。")
  }

  site_df <- readr::read_tsv(path, show_col_types = FALSE)
  missing_cols <- setdiff(required_cols, names(site_df))
  if (length(missing_cols) > 0) {
    stop("Input table is missing columns: ", paste(missing_cols, collapse = ", "))
  }
  if (nrow(site_df) == 0) stop("No selected sites found in: ", path)
  if (anyDuplicated(site_df$chr)) stop("Duplicated tRNA entries in: ", path)

  site_df <- site_df |>
    mutate(
      across(c(WT_Ratio, all_of(line_cols)), as.numeric),
      across(all_of(dep_cols), ~ .x == "TRUE"),
      family = chr |>
        stringr::str_replace("^tRNA-(Nu|Mt|Pt)-", "") |>
        stringr::str_replace("-[0-9]+[a-z]?$", "")
    )

  if (any(site_df$family == "")) {
    stop("Could not derive family label for every tRNA in: ", path)
  }

  if (any(!is.finite(site_df$WT_Ratio)) ||
      any(site_df$WT_Ratio <= 0 | site_df$WT_Ratio > 1)) {
    stop("Every selected site must have 0 < WT_Ratio <= 1.")
  }
  ratio_cols <- c("WT_Ratio", line_cols)
  if (any(!is.na(unlist(site_df[ratio_cols])) &
          (unlist(site_df[ratio_cols]) < 0 | unlist(site_df[ratio_cols]) > 1))) {
    stop("Modification ratios must lie in [0, 1] or be NA.")
  }
  for (dep_col in dep_cols) {
    ratio_col <- stringr::str_remove(dep_col, "_dependent$")
    bad <- site_df |>
      filter(!is.na(.data[[dep_col]]) & .data[[dep_col]] & is.na(.data[[ratio_col]]))
    if (nrow(bad) > 0) {
      stop("Dependent=TRUE but line ratio is NA for: ", paste(bad$chr, collapse = ", "))
    }
  }

  site_df
}

# 长表：一行一个 (位点, 样本) 的修饰比 + per-line dependent（WT 行 NA）
pivot_site_specificity_ratios <- function(site_df, genotypes) {
  line_cols <- site_specificity_line_cols(genotypes)
  dep_cols <- paste0(line_cols, "_dependent")

  ratio_long <- site_df |>
    select(chr, WT_Ratio, all_of(line_cols)) |>
    pivot_longer(cols = c("WT_Ratio", all_of(line_cols)),
                 names_to = "sample", values_to = "ratio")

  dep_long <- site_df |>
    select(chr, all_of(dep_cols)) |>
    rename_with(~ stringr::str_remove(.x, "_dependent$")) |>
    pivot_longer(cols = -chr, names_to = "sample", values_to = "is_dependent")

  ratio_long |>
    left_join(dep_long, by = c("chr", "sample")) |>
    mutate(sample = factor(sample, levels = c("WT_Ratio", line_cols)))
}

# ---------------------------------------------------------------------------
# 热图：行 = 入选 tRNA（WT 比例降序），列 = WT + 各突变体株系
# ---------------------------------------------------------------------------
run_site_specificity_heatmap <- function(input_path, genotypes, out_dir, file_prefix,
                                         plot_width = 8) {
  site_df <- prepare_site_specificity_table(input_path, genotypes = genotypes)
  line_cols <- site_specificity_line_cols(genotypes)

  # 画布与字号绑定（对齐 ../mRNA/heatmap/plot.R），行高随位点数同步放大
  base_size <- round(plot_width * 72 / 27)
  font_scale <- base_size / 16
  n_rows <- nrow(site_df)
  plot_height <- max(6, 2.4 + n_rows * 0.30 * font_scale)
  init_theme(base_size = base_size)

  sample_levels <- c("WT_Ratio", line_cols)
  sample_labels <- c("WT", paste0("*", line_cols, "*"))
  names(sample_labels) <- sample_levels

  # ggplot 的 y 因子从下到上排列：WT 比例升序，让最高值位于顶部
  row_info <- site_df |>
    select(chr, family, WT_Ratio) |>
    arrange(WT_Ratio, chr)

  plot_df <- pivot_site_specificity_ratios(site_df, genotypes) |>
    mutate(
      row_id = factor(chr, levels = row_info$chr),
      ratio_label = ifelse(
        is.na(ratio), "NA",
        ifelse(ratio < 0.10, sprintf("%.1f%%", 100 * ratio),
               sprintf("%.0f%%", 100 * ratio))
      )
    )

  p <- ggplot(plot_df, aes(x = sample, y = row_id, fill = ratio)) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_tile(
      data = filter(plot_df, is_dependent %in% TRUE),
      aes(x = sample, y = row_id),
      fill = NA, inherit.aes = FALSE,
      colour = "#D7301F", linewidth = 0.9
    ) +
    geom_text(
      aes(label = ratio_label),
      family = "Arial", fontface = "bold", colour = "black",
      size = (base_size * 0.50) / .pt, vjust = 0.5
    ) +
    scale_x_discrete(labels = sample_labels) +
    scale_y_discrete(labels = setNames(row_info$family, row_info$chr)) +
    scale_fill_gradientn(
      colours = c("white", "#E64B35CC"),
      limits = c(0, 1),
      na.value = "grey90",
      guide = "none"
    ) +
    labs(
      x = NULL, y = NULL,
      caption = "Red border: dependent in the corresponding mutant line"
    ) +
    theme(
      axis.text.x = ggtext::element_markdown(
        family = "Arial", face = "bold", colour = "black",
        size = base_size * 0.75
      ),
      axis.text.y = element_text(
        family = "Arial", face = "bold", colour = "black", size = base_size * 0.65
      ),
      axis.ticks = element_blank(),
      legend.position = "none",
      plot.caption = element_text(
        family = "Arial", face = "bold", colour = "black",
        size = base_size * 0.65, hjust = 0
      ),
      plot.caption.position = "plot",
      plot.margin = margin(10, 20, 10, 12)
    )

  output_stem <- file.path(out_dir, paste0(file_prefix, "_site_specficity_heatmap"))
  save_plot(p, output_stem, width = plot_width, height = plot_height)

  # 每基因型双株系同时 dependent 的位点数（与 final 表 AND 口径一致）
  n_dep <- vapply(genotypes, function(g) {
    sum(site_df[[paste0(g, "_1_dependent")]] &
          site_df[[paste0(g, "_2_dependent")]], na.rm = TRUE)
  }, integer(1))
  metadata <- tibble::tibble(
    n_tRNAs = nrow(site_df),
    !!!setNames(as.list(n_dep), paste0("n_dependent_", genotypes))
  )
  message("[site_specificity] heatmap: ", n_rows, " tRNAs -> ", output_stem)
  invisible(list(plot = p, metadata = metadata))
}

# ---------------------------------------------------------------------------
# 修饰比例分布：x = WT + 各株系，y = 修饰比例(%)
# ---------------------------------------------------------------------------
run_site_specificity_distribution <- function(input_path, genotypes, position_label,
                                              out_dir, file_prefix,
                                              plot_width = 7.2, plot_height = 4.8) {
  site_df <- prepare_site_specificity_table(input_path, genotypes = genotypes)
  line_cols <- site_specificity_line_cols(genotypes)
  fill_map <- site_specificity_fill_map(genotypes)

  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)

  plot_df <- pivot_site_specificity_ratios(site_df, genotypes) |>
    mutate(
      sample_group = factor(as.character(sample),
                            levels = c("WT_Ratio", line_cols),
                            labels = c("WT", line_cols)),
      ratio_pct = 100 * as.numeric(ratio)
    ) |>
    filter(!is.na(ratio_pct))

  n_na <- sum(is.na(site_df[c("WT_Ratio", line_cols)]))

  x_labels <- c("WT", paste0("<i>", line_cols, "</i>"))
  names(x_labels) <- c("WT", line_cols)

  p <- ggplot(plot_df, aes(x = sample_group, y = ratio_pct, fill = sample_group)) +
    geom_boxplot(
      width = 0.58,
      outlier.shape = NA,
      color = "black",
      linewidth = SIZE_TOKENS$base_linewidth
    ) +
    geom_point(
      position = position_jitter(width = 0.16, height = 0, seed = 1),
      shape = 21,
      size = 2.1,
      alpha = 0.95,
      color = "black",
      stroke = 0.25
    ) +
    scale_fill_manual(values = fill_map, guide = "none") +
    scale_x_discrete(labels = x_labels) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = sprintf("%s WT-modified tRNAs", position_label),
      x = NULL,
      y = paste0("Ψ", " modification level (%)")
    ) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    theme(
      axis.text.x = ggtext::element_markdown(
        family = "Arial", face = "bold", colour = "black",
        margin = margin(t = 6)
      ),
      plot.title = element_text(
        size = SIZE_TOKENS$strip_text_size,
        hjust = 0.5,
        margin = margin(b = 8)
      )
    )

  output_stem <- file.path(out_dir, paste0(file_prefix, "_site_specficity_distribution"))
  save_plot(p, output_stem, width = plot_width, height = plot_height)

  summary_df <- plot_df |>
    group_by(sample_group) |>
    summarise(
      n_values = dplyr::n(),
      median_ratio_pct = median(ratio_pct),
      .groups = "drop"
    )
  message("[site_specificity] distribution: ", nrow(site_df), " tRNAs, ",
          n_na, " NA ratio(s) dropped -> ", output_stem)
  invisible(list(plot = p, summary = summary_df))
}
