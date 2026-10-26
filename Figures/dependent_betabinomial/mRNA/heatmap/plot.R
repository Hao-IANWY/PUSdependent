#!/usr/bin/env Rscript
# plot.R — WT 修饰比例 >= 20% 的 mRNA dependent 位点跨基因型热图
#
# 输入：dependent_site_ratio_heatmap.tsv（data_prep.sh 生成）
# 行：组织内的 dependent 位点，按 WT 修饰比例由高到低排列；
# 列：WT + 7 个突变体；格值：修饰比例（两条突变株系均值），固定 0–100% 色阶；
# 格内标注百分比（<10% 保留一位小数，其余取整）；红色方框：对应突变体中判定为 dependent。
# 输出：mRNA_dependent_site_ratio_heatmap.{pdf,png,svg}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(ggtext)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  normalizePath(getwd())
}
setwd(script_dir)
shared_dir <- normalizePath(file.path(script_dir, "..", "..", "..", "shared"), mustWork = TRUE)
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

input_file <- file.path(script_dir, "dependent_site_ratio_heatmap.tsv")
if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file, "\nRun data_prep.sh first.")
}

plot_df <- read.table(
  input_file, header = TRUE, sep = "\t", quote = "", comment.char = "",
  na.strings = "NA", stringsAsFactors = FALSE, check.names = FALSE
)
required_cols <- c("tissue", "site_id", "chr", "pos", "strand", "context",
                   "WT_Ratio", "genotype", "ratio", "is_dependent")
missing_cols <- setdiff(required_cols, names(plot_df))
if (length(missing_cols) > 0) {
  stop("Input table is missing columns: ", paste(missing_cols, collapse = ", "))
}
if (nrow(plot_df) == 0) stop("No selected sites found in: ", input_file)

plot_df <- plot_df |>
  mutate(
    WT_Ratio = as.numeric(WT_Ratio),
    ratio = as.numeric(ratio),
    is_dependent = as.logical(is_dependent)
  )
if (any(!is.finite(plot_df$WT_Ratio)) ||
    any(plot_df$WT_Ratio < 0.20 | plot_df$WT_Ratio > 1)) {
  stop("Every selected site must have 0.20 <= WT_Ratio <= 1.")
}
if (any(!is.na(plot_df$ratio) & (plot_df$ratio < 0 | plot_df$ratio > 1))) {
  stop("Modification ratios must lie in [0, 1].")
}
if (anyDuplicated(plot_df[c("tissue", "site_id", "genotype")])) {
  stop("Duplicated tissue/site/genotype rows in input table.")
}

# 画布与字号绑定；行高随字号和位点数同步放大，避免长表挤压。
FIG_WIDTH <- 9
BASE_SIZE <- round(FIG_WIDTH * 72 / 27)
FONT_SCALE <- BASE_SIZE / 16
n_tissue_sites <- nrow(distinct(plot_df, tissue, site_id))
FIG_HEIGHT <- max(6, 2.4 + n_tissue_sites * 0.30 * FONT_SCALE)
init_theme(base_size = BASE_SIZE)

# 底层目录仍用 Root/leaf/Flower；图面按项目统一口径显示 Roots/Leaves/Flowers。
tissue_levels <- c("Root", "leaf", "Flower")
tissue_labels <- c("Roots", "Leaves", "Flowers")
genotype_levels <- c("WT", "trua1a", "trua1b", "trua2", "trua3", "trua4", "trua5", "trua6")
genotype_labels <- c("WT", "*trua1a*", "*trua1b*", "*trua2*",
                     "*trua3*", "*trua4*", "*trua5*", "*trua6*")
names(genotype_labels) <- genotype_levels

plot_df <- plot_df |>
  mutate(
    tissue = factor(tissue, levels = tissue_levels, labels = tissue_labels),
    genotype = factor(genotype, levels = genotype_levels),
    row_id = paste(as.character(tissue), site_id, sep = "|")
  )

# ggplot 的 y 因子从下到上排列，因此按 WT_Ratio 升序设定 levels，最高值位于顶部。
site_info <- plot_df |>
  distinct(tissue, row_id, site_id, chr, pos, WT_Ratio) |>
  arrange(tissue, WT_Ratio, site_id) |>
  mutate(site_label = paste0(chr, ":", pos))
site_labels <- setNames(site_info$site_label, site_info$row_id)
plot_df <- plot_df |>
  mutate(
    row_id = factor(row_id, levels = site_info$row_id),
    ratio_label = ifelse(
      is.na(ratio), "NA",
      ifelse(ratio < 0.10,
             sprintf("%.1f%%", 100 * ratio),
             sprintf("%.0f%%", 100 * ratio))
    )
  )

# 每个组织 × 位点必须恰好有一个 WT 值，且至少一个突变体 dependent 标记。
wt_counts <- plot_df |>
  filter(genotype == "WT") |>
  count(tissue, site_id, name = "n_wt")
if (nrow(wt_counts) != n_tissue_sites || any(wt_counts$n_wt != 1)) {
  stop("Each tissue/site must have exactly one WT ratio row.")
}
dependent_counts <- plot_df |>
  group_by(tissue, site_id) |>
  summarise(n_dependent = sum(is_dependent, na.rm = TRUE), .groups = "drop")
if (any(dependent_counts$n_dependent < 1)) {
  stop("Every selected site must be dependent in at least one mutant.")
}

p <- ggplot(plot_df, aes(x = genotype, y = row_id, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(
    data = filter(plot_df, is_dependent %in% TRUE),
    aes(x = genotype, y = row_id),
    fill = NA, inherit.aes = FALSE,
    colour = "#D7301F", linewidth = 0.6
  ) +
  geom_text(
    aes(label = ratio_label),
    family = "Arial", fontface = "bold", colour = "black",
    size = (BASE_SIZE * 0.50) / .pt, vjust = 0.5
  ) +
  facet_grid(tissue ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = genotype_labels) +
  scale_y_discrete(labels = site_labels) +
  scale_fill_gradientn(
    colours = c("white", "#E64B35CC"),
    limits = c(0, 1),
    na.value = "grey90",
    guide = "none"
  ) +
  labs(x = NULL, y = NULL, caption = "Red border: dependent in the corresponding mutant") +
  theme(
    axis.text.x = ggtext::element_markdown(
      family = "Arial", face = "bold", colour = "black",
      size = BASE_SIZE * 0.75, angle = 90, hjust = 1, vjust = 0.5
    ),
    axis.text.y = element_text(
      family = "Arial", face = "bold", colour = "black", size = BASE_SIZE * 0.65
    ),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      family = "Arial", face = "bold", colour = "black", size = BASE_SIZE * 0.80
    ),
    panel.spacing.y = grid::unit(0.5, "lines"),
    legend.position = "none",
    plot.caption = element_text(
      family = "Arial", face = "bold", colour = "black",
      size = BASE_SIZE * 0.65, hjust = 0
    ),
    plot.caption.position = "plot",
    plot.margin = margin(10, 20, 10, 12)
  )

output_stem <- file.path(script_dir, "mRNA_dependent_site_ratio_heatmap")
save_plot(p, output_stem, width = FIG_WIDTH, height = FIG_HEIGHT)
message("[plot] Done: ", n_tissue_sites,
        " tissue-site rows; PDF + 600-DPI PNG + SVG written to ", script_dir)
