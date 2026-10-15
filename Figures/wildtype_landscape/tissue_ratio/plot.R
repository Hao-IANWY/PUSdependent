#!/usr/bin/env Rscript
# WT positive 位点修饰率的三组织分布比较(mRNA / rRNA / tRNA)。
# 柱 = 位点均值,误差棒 = SEM,散点 = 单位点修饰率,括号 = Welch t 检验。
# 参考 archive/20260902/paper/WT_sites_overview/tissue_ratio,
# 数据源与样式对齐当前仓库(polyARNA + tRNA 双流程,Figures/shared 样式体系)。

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]

  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }

  normalizePath(getwd())
}

script_dir <- get_script_dir()
repo_dir <- normalizePath(file.path(script_dir, "..", "..", ".."))
source(file.path(repo_dir, "Figures/shared/plot_style.R"))
source(file.path(repo_dir, "Figures/shared/plot_tokens.R"))

init_theme(base_size = 16, border = TRUE, sizes = SIZE_TOKENS)

tissue_levels <- c("Root", "Leaf", "Flower")
rna_levels <- c("mRNA", "rRNA", "tRNA")

# polyARNA 样本名连字符(WT-leaf),tRNA 样本名下划线(WT_Leaf)
positive_specs <- tribble(
  ~rna_type, ~tissue, ~positive_file, ~fisher_file, ~chr_regex,
  "mRNA", "Root",   "polyARNA/STAT/fisher/mRNA/WT-Root_sites.positive.tsv",   "polyARNA/STAT/fisher/mRNA/WT-Root_sites.Fisher.tsv",   NA_character_,
  "mRNA", "Leaf",   "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv",   "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.Fisher.tsv",   NA_character_,
  "mRNA", "Flower", "polyARNA/STAT/fisher/mRNA/WT-Flower_sites.positive.tsv", "polyARNA/STAT/fisher/mRNA/WT-Flower_sites.Fisher.tsv", NA_character_,
  "rRNA", "Root",   "polyARNA/STAT/fisher/rRNA/WT-Root_sites.positive.tsv",   "polyARNA/STAT/fisher/rRNA/WT-Root_sites.Fisher.tsv",   NA_character_,
  "rRNA", "Leaf",   "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv",   "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.Fisher.tsv",   NA_character_,
  "rRNA", "Flower", "polyARNA/STAT/fisher/rRNA/WT-Flower_sites.positive.tsv", "polyARNA/STAT/fisher/rRNA/WT-Flower_sites.Fisher.tsv", NA_character_,
  "tRNA", "Root",   "tRNA/STAT/fisher/ncRNA/WT_Root_sites.positive.tsv",   "tRNA/STAT/fisher/ncRNA/WT_Root_sites.Fisher.tsv",   "^tRNA-",
  "tRNA", "Leaf",   "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv",   "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.Fisher.tsv",   "^tRNA-",
  "tRNA", "Flower", "tRNA/STAT/fisher/ncRNA/WT_Flower_sites.positive.tsv", "tRNA/STAT/fisher/ncRNA/WT_Flower_sites.Fisher.tsv", "^tRNA-"
) %>%
  mutate(
    positive_file = file.path(repo_dir, positive_file),
    fisher_file = file.path(repo_dir, fisher_file)
  )

require_columns <- function(tbl, cols, file_path) {
  missing_cols <- setdiff(cols, colnames(tbl))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }
}

read_site_data <- function(rna_type, tissue, positive_file, fisher_file, chr_regex) {
  positive <- readr::read_tsv(positive_file, show_col_types = FALSE)
  fisher <- readr::read_tsv(fisher_file, show_col_types = FALSE)

  require_columns(positive, c("chr", "pos", "strand", "ratio", "depth"), positive_file)
  require_columns(fisher, c("chr", "pos", "strand", "Treat_avg_depth", "Treat_avg_gap"), fisher_file)

  positive <- positive %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand)
    )

  fisher <- fisher %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      Treat_avg_depth = as.numeric(Treat_avg_depth),
      Treat_avg_gap = as.numeric(Treat_avg_gap)
    )

  # tRNA 只统计核编码 tRNA,排除 sn/snoRNA 等同参考映射记录
  if (!is.na(chr_regex) && nzchar(chr_regex)) {
    positive <- positive %>% filter(str_detect(chr, chr_regex))
    fisher <- fisher %>% filter(str_detect(chr, chr_regex))
  }

  positive %>%
    left_join(
      fisher %>% select(chr, pos, strand, Treat_avg_depth, Treat_avg_gap),
      by = c("chr", "pos", "strand")
    ) %>%
    transmute(
      rna_type = factor(rna_type, levels = rna_levels),
      tissue = factor(tissue, levels = tissue_levels),
      chr,
      pos,
      strand,
      site_id = paste(chr, pos, strand, sep = ":"),
      ratio = as.numeric(ratio),
      ratio_pct = ratio * 100,
      depth = as.numeric(depth),
      treat_depth = dplyr::coalesce(Treat_avg_depth, depth),
      treat_gap = dplyr::coalesce(Treat_avg_gap, ratio * depth)
    ) %>%
    filter(
      is.finite(ratio_pct),
      ratio_pct >= 0,
      ratio_pct <= 100,
      is.finite(treat_depth),
      treat_depth > 0,
      is.finite(treat_gap),
      treat_gap >= 0
    )
}

site_tbl <- pmap_dfr(positive_specs, read_site_data) %>%
  mutate(
    rna_type = factor(as.character(rna_type), levels = rna_levels),
    tissue = factor(as.character(tissue), levels = tissue_levels)
  )

if (nrow(site_tbl) == 0) {
  stop("No WT positive sites were loaded.")
}

count_tbl <- site_tbl %>%
  mutate(
    modified_reads = pmax(0, round(treat_gap)),
    total_reads = pmax(0, round(treat_depth)),
    unmodified_reads = pmax(0, total_reads - modified_reads)
  ) %>%
  group_by(rna_type, tissue) %>%
  summarise(
    n_sites = n(),
    median_ratio_pct = median(ratio_pct, na.rm = TRUE),
    mean_ratio_pct = mean(ratio_pct, na.rm = TRUE),
    modified_reads = sum(modified_reads, na.rm = TRUE),
    unmodified_reads = sum(unmodified_reads, na.rm = TRUE),
    total_reads = modified_reads + unmodified_reads,
    aggregate_ratio_pct = 100 * modified_reads / total_reads,
    .groups = "drop"
  )

# ---- 组间 Welch t 检验(位点级 ratio_pct,不校正多重检验) ----
pair_specs <- tribble(
  ~tissue_a, ~tissue_b,
  "Root", "Leaf",
  "Root", "Flower",
  "Leaf", "Flower"
)

pvalue_to_star <- function(p) {
  case_when(
    is.na(p) ~ "n.s.",
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

safe_t_test <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  fallback <- tibble(
    p_value = NA_real_,
    statistic = NA_real_,
    parameter = NA_real_,
    mean_a = mean(x, na.rm = TRUE),
    mean_b = mean(y, na.rm = TRUE),
    median_a = median(x, na.rm = TRUE),
    median_b = median(y, na.rm = TRUE)
  )

  if (length(x) < 2 || length(y) < 2 || length(unique(c(x, y))) < 2) {
    return(fallback)
  }

  test <- tryCatch(
    stats::t.test(x, y, var.equal = FALSE),
    error = function(e) NULL
  )

  if (is.null(test)) {
    return(fallback)
  }

  tibble(
    p_value = unname(test$p.value),
    statistic = unname(test$statistic),
    parameter = unname(test$parameter),
    mean_a = mean(x, na.rm = TRUE),
    mean_b = mean(y, na.rm = TRUE),
    median_a = median(x, na.rm = TRUE),
    median_b = median(y, na.rm = TRUE)
  )
}

ttest_tbl <- tidyr::crossing(rna_type = factor(rna_levels, levels = rna_levels), pair_specs) %>%
  rowwise() %>%
  mutate(
    test_result = list(safe_t_test(
      site_tbl %>%
        filter(rna_type == .env$rna_type, tissue == .env$tissue_a) %>%
        pull(ratio_pct),
      site_tbl %>%
        filter(rna_type == .env$rna_type, tissue == .env$tissue_b) %>%
        pull(ratio_pct)
    ))
  ) %>%
  ungroup() %>%
  unnest(test_result) %>%
  left_join(
    count_tbl %>%
      select(rna_type, tissue_a = tissue, n_sites_a = n_sites),
    by = c("rna_type", "tissue_a")
  ) %>%
  left_join(
    count_tbl %>%
      select(rna_type, tissue_b = tissue, n_sites_b = n_sites),
    by = c("rna_type", "tissue_b")
  ) %>%
  mutate(
    significance = pvalue_to_star(p_value),
    comparison = paste(tissue_a, tissue_b, sep = " vs "),
    delta_mean = mean_a - mean_b,
    delta_median = median_a - median_b
  ) %>%
  select(
    rna_type, comparison, tissue_a, tissue_b,
    n_sites_a, n_sites_b,
    mean_a, mean_b, delta_mean,
    median_a, median_b, delta_median,
    statistic, parameter,
    p_value, significance
  )

# ---- 绘图数据:组内三组织横向错位,显著性括号架在组上方 ----
plot_tbl <- site_tbl %>%
  mutate(
    rna_index = as.integer(rna_type),
    tissue_index = as.integer(tissue),
    x = rna_index + (tissue_index - 2) * 0.26,
    x_group = factor(paste(rna_type, tissue, sep = "_"),
                     levels = as.vector(outer(rna_levels, tissue_levels, paste, sep = "_")))
  )

bar_tbl <- plot_tbl %>%
  group_by(rna_type, tissue, x_group, x) %>%
  summarise(
    mean_ratio_pct = mean(ratio_pct, na.rm = TRUE),
    sd_ratio_pct = stats::sd(ratio_pct, na.rm = TRUE),
    n_sites = n(),
    sem_ratio_pct = sd_ratio_pct / sqrt(n_sites),
    .groups = "drop"
  )

annotation_base <- plot_tbl %>%
  group_by(rna_type) %>%
  summarise(y_max = max(ratio_pct, na.rm = TRUE), .groups = "drop")

annotation_tbl <- ttest_tbl %>%
  mutate(
    rna_type = factor(as.character(rna_type), levels = rna_levels),
    pair_index = match(comparison, c("Root vs Leaf", "Root vs Flower", "Leaf vs Flower")),
    x_base = as.integer(rna_type),
    x_start = x_base + (match(tissue_a, tissue_levels) - 2) * 0.26,
    x_end = x_base + (match(tissue_b, tissue_levels) - 2) * 0.26
  ) %>%
  left_join(annotation_base, by = "rna_type") %>%
  mutate(
    y = pmax(106, y_max + 2) + pair_index * 10,
    y_text = y + 2.5,
    label = significance
  )

# Y 轴上限罩住最高刻度:把括号文字层也算进被画的 y 值
y_ax <- y_axis_cover(c(plot_tbl$ratio_pct, annotation_tbl$y, annotation_tbl$y_text))

tissue_colors <- setNames(PAL_TISSUE, tissue_levels)

ratio_plot <- ggplot() +
  geom_col(
    data = bar_tbl,
    aes(x = x, y = mean_ratio_pct, fill = tissue, group = x_group),
    width = 0.20,
    alpha = 0.58,
    color = "black",
    linewidth = SIZE_TOKENS$base_linewidth
  ) +
  geom_errorbar(
    data = bar_tbl,
    aes(
      x = x,
      ymin = pmax(0, mean_ratio_pct - sem_ratio_pct),
      ymax = pmin(100, mean_ratio_pct + sem_ratio_pct),
      group = x_group
    ),
    width = 0.08,
    color = "black",
    linewidth = SIZE_TOKENS$base_linewidth
  ) +
  geom_point(
    data = plot_tbl,
    aes(x = x, y = ratio_pct, fill = tissue, group = x_group),
    position = position_jitter(width = 0.045, height = 0, seed = 20261015),
    shape = 21,
    size = 1.7,
    stroke = 0.35,
    alpha = 0.72,
    color = "black",
    show.legend = FALSE
  ) +
  geom_segment(
    data = annotation_tbl,
    aes(x = x_start, xend = x_end, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = SIZE_TOKENS$base_linewidth,
    color = "black"
  ) +
  geom_segment(
    data = annotation_tbl,
    aes(x = x_start, xend = x_start, y = y - 2.4, yend = y),
    inherit.aes = FALSE,
    linewidth = SIZE_TOKENS$base_linewidth,
    color = "black"
  ) +
  geom_segment(
    data = annotation_tbl,
    aes(x = x_end, xend = x_end, y = y - 2.4, yend = y),
    inherit.aes = FALSE,
    linewidth = SIZE_TOKENS$base_linewidth,
    color = "black"
  ) +
  geom_text(
    data = annotation_tbl,
    aes(x = (x_start + x_end) / 2, y = y_text, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 5.2,
    color = "black"
  ) +
  scale_x_continuous(
    breaks = seq_along(rna_levels),
    labels = rna_levels,
    limits = c(0.45, length(rna_levels) + 0.55),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  scale_y_continuous(
    limits = c(0, y_ax$top),
    breaks = y_ax$breaks,
    expand = y_ax$expand
  ) +
  # 显示层用复数(仓库口径),数据层 tissue 保持单数
  scale_fill_manual(
    values = tissue_colors,
    name = NULL,
    breaks = tissue_levels,
    labels = c("Roots", "Leaves", "Flowers"),
    drop = FALSE
  ) +
  labs(x = NULL, y = expression(Psi~"modification ratio (%)")) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key.width = grid::unit(18, "pt"),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = SIZE_TOKENS$axis_text_size),
    axis.title.y = element_text(size = SIZE_TOKENS$axis_title_size, margin = margin(r = 8)),
    plot.margin = margin(8, 16, 8, 12)
  )

outputs <- list(
  sites = file.path(script_dir, "tissue_ratio_site_ratios.tsv"),
  summary = file.path(script_dir, "tissue_ratio_summary.tsv"),
  ttest = file.path(script_dir, "tissue_ratio_ttest_tests.tsv"),
  plot = file.path(script_dir, "tissue_ratio_distribution")
)

readr::write_tsv(
  site_tbl %>%
    arrange(rna_type, tissue, chr, pos, strand) %>%
    select(rna_type, tissue, chr, pos, strand, site_id, ratio, ratio_pct, depth, treat_depth, treat_gap),
  outputs$sites
)
readr::write_tsv(count_tbl %>% arrange(rna_type, tissue), outputs$summary)
readr::write_tsv(ttest_tbl %>% arrange(rna_type, comparison), outputs$ttest)

save_plot(ratio_plot, outputs$plot, width = 7.2, height = 5.2)

message("Wrote ", outputs$sites)
message("Wrote ", outputs$summary)
message("Wrote ", outputs$ttest)
