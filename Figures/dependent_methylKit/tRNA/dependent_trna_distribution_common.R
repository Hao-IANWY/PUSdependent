#!/usr/bin/env Rscript
# dependent_trna_distribution_common.R — 单个基因型的 dependent 位点修饰比例分布箱线图
#
# 移植自 archive/20260902/paper/tRNA_dependent/dependent_trna_distribution_common.R，
# 样式改走 Figures/shared（plot_style.R + plot_tokens.R），导出改走 save_plot()。
#
# 每个 panel 一个 (compartment, fix_pos) 组合：x = 组织，y = Ψ 修饰比例 (%)，
# 分成 WT / line1 / line2 三组（line1/line2 = 该基因型两个株系各自 3 个 Treat 重复的
# Mut_Ratio 均值）。数据来自 data_prep.sh 生成的 <genotype>_all_tissues.dependent.tsv。

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(ggtext)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

distribution_colors <- c(
  WT = "#9E9E9ECC",
  line1 = "#F4A3AECC",
  line2 = "#E26D7FCC"
)

distribution_tissue_levels <- c("Roots", "Leaves", "Flowers")
distribution_tissue_source_levels <- c("Root", "leaf", "Flower")

build_distribution_theme <- function() {
  theme_base(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS) +
    ggplot2::theme(
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10))
    )
}

# 图例标签：WT 正体，两个株系名用斜体（ggtext markdown 渲染为 Arial Bold Italic）
build_distribution_legend_labels <- function(genotype) {
  c(
    "WT",
    paste0("<i>", genotype, "_1</i>"),
    paste0("<i>", genotype, "_2</i>")
  )
}

setup_distribution_style <- function() {
  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)
}

normalize_distribution_specs <- function(plot_specs) {
  plot_specs |>
    mutate(
      compartment = as.character(compartment),
      label = as.character(label),
      fix_pos_values = purrr::map(fix_pos_values, as.integer)
    )
}

prepare_single_panel_distribution_data <- function(site_df, genotype, plot_specs) {
  line_cols <- c(paste0(genotype, "_1"), paste0(genotype, "_2"))
  required_cols <- c("chr", "pos", "fix_pos", "WT_Ratio", line_cols, "Tissue")
  missing_cols <- setdiff(required_cols, names(site_df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in dependent-site table: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_specs <- normalize_distribution_specs(plot_specs)

  base_df <- site_df |>
    mutate(
      Tissue = factor(
        as.character(Tissue),
        levels = distribution_tissue_source_levels,
        labels = distribution_tissue_levels
      ),
      chr = as.character(chr),
      pos = as.integer(pos),
      fix_pos = as.integer(fix_pos),
      compartment = dplyr::case_when(
        stringr::str_starts(chr, "tRNA-Nu-") ~ "Cy",
        stringr::str_starts(chr, "tRNA-Mt-") ~ "Mt",
        stringr::str_starts(chr, "tRNA-Pt-") ~ "Pt",
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(Tissue), !is.na(compartment))

  filtered_df <- purrr::map_dfr(seq_len(nrow(plot_specs)), function(i) {
    spec <- plot_specs[i, ]
    base_df |>
      filter(
        compartment == spec$compartment[[1]],
        fix_pos %in% spec$fix_pos_values[[1]]
      ) |>
      mutate(position_group = spec$label[[1]])
  })

  if (nrow(filtered_df) == 0) {
    stop("No dependent sites matched the requested compartment/position groups for ", genotype, ".")
  }

  filtered_df |>
    dplyr::select(Tissue, chr, pos, fix_pos, compartment, position_group,
                  WT_Ratio, dplyr::all_of(line_cols)) |>
    tidyr::pivot_longer(
      cols = c("WT_Ratio", dplyr::all_of(line_cols)),
      names_to = "sample",
      values_to = "ratio"
    ) |>
    mutate(
      sample_group = dplyr::case_when(
        sample == "WT_Ratio" ~ "WT",
        sample == line_cols[1] ~ "line1",
        sample == line_cols[2] ~ "line2",
        TRUE ~ sample
      ),
      ratio = as.numeric(ratio),
      ratio_pct = ratio * 100,
      position_group = factor(position_group, levels = plot_specs$label),
      sample_group = factor(sample_group, levels = c("WT", "line1", "line2"))
    ) |>
    filter(
      !is.na(position_group),
      !is.na(ratio_pct),
      ratio_pct >= 0,
      ratio_pct <= 100
    ) |>
    arrange(position_group, Tissue, sample_group, chr, pos)
}

build_single_panel_distribution_plot <- function(plot_df, genotype) {
  dodge_width <- 0.72
  if (nrow(plot_df) == 0) {
    stop("No rows remained for the single-panel position-distribution plot.")
  }

  position_labels <- unique(as.character(plot_df$position_group))
  if (length(position_labels) != 1L) {
    stop("Each single-panel distribution plot must contain exactly one position group.")
  }

  ggplot(plot_df, aes(x = Tissue, y = ratio_pct, fill = sample_group)) +
    geom_boxplot(
      position = position_dodge(width = dodge_width),
      width = 0.58,
      outlier.shape = NA,
      color = "black",
      linewidth = SIZE_TOKENS$base_linewidth
    ) +
    geom_point(
      aes(fill = sample_group),
      position = position_jitterdodge(
        jitter.width = 0.14,
        jitter.height = 0,
        dodge.width = dodge_width,
        seed = 1
      ),
      shape = 21,
      size = 2.1,
      alpha = 0.95,
      color = "black",
      stroke = 0.25
    ) +
    scale_fill_manual(
      values = distribution_colors,
      breaks = c("WT", "line1", "line2"),
      labels = build_distribution_legend_labels(genotype)
    ) +
    scale_x_discrete(drop = FALSE, limits = distribution_tissue_levels) +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = position_labels[[1]],
      x = NULL,
      y = paste0("Ψ", " modification level (%)"),
      fill = NULL
    ) +
    coord_cartesian(ylim = c(0, 100), clip = "off") +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.5, 0.97),
      legend.justification = c(0.5, 1),
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.background = element_blank(),
      legend.background = element_blank(),
      legend.text = ggtext::element_markdown(size = 15),
      legend.key.width = grid::unit(5, "mm"),
      legend.key.height = grid::unit(3.8, "mm"),
      legend.spacing.x = grid::unit(2, "mm"),
      legend.margin = margin(1, 0, 1, 0),
      axis.text.x = element_text(margin = margin(t = 6)),
      plot.title = element_text(
        size = SIZE_TOKENS$strip_text_size,
        hjust = 0.5,
        margin = margin(b = 8)
      )
    )
}

remove_obsolete_distribution_outputs <- function(out_dir, genotype) {
  obsolete_files <- list.files(
    out_dir,
    pattern = paste0("^", genotype, ".*dependent_position_distribution\\.(pdf|png|svg)$"),
    full.names = TRUE
  )

  invisible(file.remove(obsolete_files[file.exists(obsolete_files)]))
}

save_single_panel_distribution_plots <- function(
  plot_df,
  out_dir,
  genotype,
  plot_specs,
  save_outputs = TRUE,
  plot_width = 6.2,
  plot_height = 4.8
) {
  plot_specs <- normalize_distribution_specs(plot_specs) |>
    mutate(output_suffix = as.character(output_suffix))
  plot_data_tsv <- file.path(out_dir, paste0(genotype, "_dependent_position_distribution.tsv"))
  summary_tsv <- file.path(out_dir, paste0(genotype, "_dependent_position_distribution_summary.tsv"))

  summary_df <- plot_df |>
    group_by(position_group, Tissue, sample_group) |>
    summarise(
      n_values = dplyr::n(),
      median_ratio_pct = median(ratio_pct, na.rm = TRUE),
      mean_ratio_pct = mean(ratio_pct, na.rm = TRUE),
      .groups = "drop"
    )

  output_tbl <- plot_specs |>
    transmute(
      position_group = factor(label, levels = plot_specs$label),
      plot_pdf = NA_character_,
      plot_png = NA_character_,
      plot_svg = NA_character_
    )

  if (save_outputs) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    remove_obsolete_distribution_outputs(out_dir, genotype = genotype)
    readr::write_tsv(plot_df, plot_data_tsv)
    readr::write_tsv(summary_df, summary_tsv)

    output_tbl <- purrr::map_dfr(seq_len(nrow(plot_specs)), function(i) {
      spec <- plot_specs[i, ]
      position_label <- spec$label[[1]]
      group_df <- plot_df |>
        filter(as.character(position_group) == position_label) |>
        droplevels()
      p <- build_single_panel_distribution_plot(group_df, genotype = genotype)
      output_stem <- file.path(
        out_dir,
        paste0(genotype, "_", spec$output_suffix[[1]], "_dependent_position_distribution")
      )
      save_plot(
        plot = p,
        filename = output_stem,
        width = plot_width,
        height = plot_height,
        formats = c("pdf", "png", "svg"),
        dpi = export_dpi
      )

      tibble::tibble(
        position_group = factor(position_label, levels = plot_specs$label),
        plot_pdf = paste0(output_stem, ".pdf"),
        plot_png = paste0(output_stem, ".png"),
        plot_svg = paste0(output_stem, ".svg")
      )
    })
  }

  metadata <- summary_df |>
    group_by(position_group) |>
    summarise(
      n_values = sum(n_values),
      n_tissues = dplyr::n_distinct(Tissue),
      .groups = "drop"
    ) |>
    left_join(output_tbl, by = "position_group")

  list(
    plot_data = plot_df,
    summary = summary_df,
    metadata = metadata,
    plot_data_tsv = plot_data_tsv,
    summary_tsv = summary_tsv
  )
}

run_single_panel_dependent_position_distributions <- function(
  genotype,
  input_path,
  out_dir,
  plot_specs,
  save_outputs = TRUE,
  plot_width = 6.2,
  plot_height = 4.8
) {
  if (!file.exists(input_path)) {
    stop(
      "Missing input table: ", input_path,
      "\n先跑 data_prep.sh 生成。"
    )
  }

  setup_distribution_style()
  plot_df <- prepare_single_panel_distribution_data(
    site_df = readr::read_tsv(input_path, show_col_types = FALSE),
    genotype = genotype,
    plot_specs = plot_specs
  )

  save_single_panel_distribution_plots(
    plot_df = plot_df,
    out_dir = out_dir,
    genotype = genotype,
    plot_specs = plot_specs,
    save_outputs = save_outputs,
    plot_width = plot_width,
    plot_height = plot_height
  )
}
