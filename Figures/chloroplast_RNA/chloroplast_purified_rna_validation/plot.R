#!/usr/bin/env Rscript
# chloroplast_purified_rna_validation — 叶绿体纯化 RNA 建库的 Ψ 位点数量与重叠验证
#
# 移植自 archive/20260902/paper/chloroplast_RNA/chloroplast_purified_rna_validation/plot.R。
# 数据源换成当前流程：
#   Chloroplast/STAT/fisher/ncRNA/Chloroplast_sites.positive.tsv  (rRNA Ch23S / tRNA ^tRNA-Pt-)
#   Chloroplast/STAT/fisher/mRNA/Chloroplast_sites.positive.tsv   (mRNA ^ATCG)
# 对照的 WT leaf 位点直接复用 chloroplast_site_overview 的产出表。
#
# rRNA / mRNA 与 poly(A)+ 数据集比较，tRNA 与 small RNA 数据集比较；
# overlap 以 chr:pos:strand 为位点键。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
  library(tibble)
})

rna_levels <- c("rRNA", "tRNA", "mRNA")

count_fill <- c(
  "rRNA" = "#8DD3C7",
  "tRNA" = "#B3DE69",
  "mRNA" = "#FDB462"
)

venn_fill <- c(
  "Leaf" = "#7FB6D6",
  "Chloroplast RNA" = "#E9A35C"
)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]
  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }
  normalizePath(getwd())
}

read_chloroplast_sites <- function(file_path, rna_type, chr_pattern) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "ratio", "depth")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      ratio = as.numeric(ratio),
      depth = as.numeric(depth)
    ) %>%
    filter(grepl(chr_pattern, chr)) %>%
    transmute(
      rna_type = factor(rna_type, levels = rna_levels),
      chr,
      pos,
      strand,
      ratio,
      depth,
      ratio_pct = ratio * 100,
      site_id = paste(chr, pos, strand, sep = ":")
    ) %>%
    filter(is.finite(ratio_pct), ratio_pct >= 0, ratio_pct <= 100)
}

circle_points <- function(center_x, center_y, radius, n = 240) {
  theta <- seq(0, 2 * pi, length.out = n)
  tibble(
    x = center_x + radius * cos(theta),
    y = center_y + radius * sin(theta)
  )
}

run_analysis <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)

  input_specs <- list(
    list(
      file = file.path(repo_root, "Chloroplast/STAT/fisher/ncRNA/Chloroplast_sites.positive.tsv"),
      rna_type = "rRNA",
      chr_pattern = "^Ch23S$"
    ),
    list(
      file = file.path(repo_root, "Chloroplast/STAT/fisher/ncRNA/Chloroplast_sites.positive.tsv"),
      rna_type = "tRNA",
      chr_pattern = "^tRNA-Pt-"
    ),
    list(
      file = file.path(repo_root, "Chloroplast/STAT/fisher/mRNA/Chloroplast_sites.positive.tsv"),
      rna_type = "mRNA",
      chr_pattern = "^ATCG"
    )
  )

  leaf_sites_file <- file.path(
    script_dir, "..", "chloroplast_site_overview",
    "chloroplast_site_overview_density_sites.tsv"
  )

  missing_inputs <- c(
    vapply(input_specs, function(x) !file.exists(x$file), logical(1)),
    !file.exists(leaf_sites_file)
  )
  if (any(missing_inputs)) {
    candidate_files <- c(vapply(input_specs, `[[`, character(1), "file"), leaf_sites_file)
    stop("Missing required inputs:\n", paste(candidate_files[missing_inputs], collapse = "\n"))
  }

  site_tbl <- bind_rows(lapply(input_specs, function(spec) {
    read_chloroplast_sites(spec$file, spec$rna_type, spec$chr_pattern)
  })) %>%
    arrange(rna_type, chr, pos, strand)

  if (nrow(site_tbl) == 0) {
    stop("No chloroplast RNA sites were found after filtering.")
  }

  leaf_sites <- readr::read_tsv(leaf_sites_file, show_col_types = FALSE) %>%
    mutate(
      rna_type = factor(as.character(rna_type), levels = rna_levels),
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      site_id = paste(chr, pos, strand, sep = ":")
    ) %>%
    distinct(rna_type, site_id)

  count_tbl <- site_tbl %>%
    count(rna_type, name = "n_sites") %>%
    complete(rna_type = factor(rna_levels, levels = rna_levels), fill = list(n_sites = 0)) %>%
    mutate(label = scales::comma(n_sites))

  overlap_summary <- bind_rows(lapply(rna_levels, function(type_name) {
    purified_ids <- site_tbl %>%
      filter(rna_type == type_name) %>%
      pull(site_id) %>%
      unique()
    leaf_ids <- leaf_sites %>%
      filter(rna_type == type_name) %>%
      pull(site_id) %>%
      unique()

    tibble(
      rna_type = factor(type_name, levels = rna_levels),
      purified_sites = length(purified_ids),
      leaf_sites = length(leaf_ids),
      overlap_sites = length(intersect(purified_ids, leaf_ids)),
      purified_only_sites = length(setdiff(purified_ids, leaf_ids)),
      leaf_only_sites = length(setdiff(leaf_ids, purified_ids))
    )
  })) %>%
    mutate(
      purified_overlap_pct = if_else(purified_sites > 0, overlap_sites / purified_sites * 100, NA_real_),
      leaf_overlap_pct = if_else(leaf_sites > 0, overlap_sites / leaf_sites * 100, NA_real_)
    )

  readr::write_tsv(site_tbl, file.path(script_dir, "chloroplast_purified_rna_validation_sites.tsv"))
  readr::write_tsv(count_tbl, file.path(script_dir, "chloroplast_purified_rna_validation_count_summary.tsv"))
  readr::write_tsv(overlap_summary, file.path(script_dir, "chloroplast_purified_rna_validation_overlap_summary.tsv"))

  # ---------------------------------------------------------------------------
  # 1. 位点数量 barplot
  # ---------------------------------------------------------------------------
  # 顶端刻度由 pretty() 保证罩住最高柱；数据外再留 1.12 倍作为刻度上界，
  # 只为柱顶数字留出文字空间。expand 仍按真实数据最大值算，面板顶正好落在最高刻度。
  n_max <- max(count_tbl$n_sites)
  y_breaks <- pretty(c(0, n_max * 1.12), n = 4)
  y_top <- max(y_breaks)
  y_expand <- expansion(mult = c(0.02, (y_top - n_max) / n_max))

  barplot <- ggplot(count_tbl, aes(x = rna_type, y = n_sites, fill = rna_type)) +
    geom_col(width = 0.68, color = "black", linewidth = SIZE_TOKENS$base_linewidth) +
    geom_text(
      aes(label = label),
      vjust = -0.45,
      family = "Arial",
      fontface = "bold",
      size = SIZE_TOKENS$annotation_text_size / ggplot2::.pt
    ) +
    scale_fill_manual(values = count_fill, drop = FALSE) +
    scale_y_continuous(
      name = expression("Number of " * Psi * " sites"),
      breaks = y_breaks,
      expand = y_expand
    ) +
    scale_x_discrete(name = NULL, drop = FALSE) +
    labs(title = "Chloroplast RNA") +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = SIZE_TOKENS$plot_title_size)
    )

  # ---------------------------------------------------------------------------
  # 2. 三个 2-set Venn（手绘圆，避免 VennDiagram 空集合报错）
  # ---------------------------------------------------------------------------
  venn_circle_tbl <- bind_rows(lapply(rna_levels, function(type_name) {
    bind_rows(
      circle_points(-0.45, 0, 0.85) %>%
        mutate(rna_type = type_name, set_name = "Leaf"),
      circle_points(0.45, 0, 0.85) %>%
        mutate(rna_type = type_name, set_name = "Chloroplast RNA")
    )
  })) %>%
    mutate(
      rna_type = factor(rna_type, levels = rna_levels),
      set_name = factor(set_name, levels = c("Leaf", "Chloroplast RNA"))
    )

  venn_count_tbl <- overlap_summary %>%
    transmute(rna_type, leaf_only_sites, overlap_sites, purified_only_sites) %>%
    pivot_longer(
      cols = c(leaf_only_sites, overlap_sites, purified_only_sites),
      names_to = "region",
      values_to = "n_sites"
    ) %>%
    mutate(
      x = case_when(
        region == "leaf_only_sites" ~ -0.68,
        region == "overlap_sites" ~ 0,
        region == "purified_only_sites" ~ 0.68
      ),
      y = 0,
      label = as.character(n_sites)
    )

  venn_set_label_tbl <- bind_rows(
    tibble(
      rna_type = factor(rna_levels, levels = rna_levels),
      set_name = factor("Leaf", levels = c("Leaf", "Chloroplast RNA")),
      label = c("poly(A)+ RNA", "small RNA", "poly(A)+ RNA")
    ),
    tibble(
      rna_type = factor(rna_levels, levels = rna_levels),
      set_name = factor("Chloroplast RNA", levels = c("Leaf", "Chloroplast RNA")),
      label = "Chloroplast RNA"
    )
  ) %>%
    mutate(
      x = if_else(set_name == "Leaf", -0.92, 0.92),
      y = -1.05
    )

  overlap_plot <- ggplot() +
    geom_polygon(
      data = venn_circle_tbl,
      aes(x = x, y = y, group = interaction(rna_type, set_name), fill = set_name),
      color = "black",
      linewidth = SIZE_TOKENS$base_linewidth,
      alpha = 0.58
    ) +
    geom_text(
      data = venn_count_tbl,
      aes(x = x, y = y, label = label),
      family = "Arial",
      fontface = "bold",
      size = SIZE_TOKENS$legend_text_size / ggplot2::.pt,
      color = "black"
    ) +
    geom_text(
      data = venn_set_label_tbl,
      aes(x = x, y = y, label = label),
      family = "Arial",
      fontface = "bold",
      size = (SIZE_TOKENS$annotation_text_size - 2) / ggplot2::.pt,
      color = "black"
    ) +
    facet_wrap(~rna_type, nrow = 1) +
    scale_fill_manual(values = venn_fill, drop = FALSE, name = NULL) +
    labs(title = "Overlap with poly(A)+ and small RNA datasets") +
    coord_fixed(xlim = c(-1.65, 1.65), ylim = c(-1.25, 1.15), clip = "off") +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = SIZE_TOKENS$axis_title_size, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = SIZE_TOKENS$plot_title_size)
    )

  save_plot(barplot, file.path(script_dir, "chloroplast_purified_rna_validation_barplot"), width = 8, height = 6.6)
  save_plot(overlap_plot, file.path(script_dir, "chloroplast_purified_rna_validation_overlap"), width = 12, height = 4.8)
}

if (sys.nframe() == 0) {
  run_analysis()
}
