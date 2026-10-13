#!/usr/bin/env Rscript
# ============================================================================
# plot.R — mRNA TRUA-dependent 位点数堆叠柱状图（XY offset 分离轴）
#
# 输入：dependent_site_count.tsv  — genotype | tissue | count（data_prep.sh 生成）
# 输出：dependent_site_count_stacked_barplot.pdf / .png / .svg
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(dplyr)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

# ---- 唯一数据来源：data_prep.sh 生成的 TSV ----
df <- read.table("dependent_site_count.tsv", header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)

df[["genotype"]] <- factor(
  df[["genotype"]],
  levels = c("trua1a", "trua1b", "trua2", "trua3", "trua4", "trua5", "trua6")
)
df[["tissue"]] <- df[["tissue"]] |>
  recode(Root = "Roots", leaf = "Leaves", Flower = "Flowers") |>
  factor(levels = c("Roots", "Leaves", "Flowers"))

plot_df <- df |>
  complete(genotype, tissue, fill = list(count = 0))

fill_colors <- setNames(PAL_TISSUE, c("Roots", "Leaves", "Flowers"))

# offset 分离轴（去闭合边框）
init_theme_offset(base_size = 24)

p <- ggplot(plot_df, aes(x = .data[["genotype"]], y = .data[["count"]],
                          fill = .data[["tissue"]])) +
  geom_col(width = 0.72, linewidth = 0.6, color = "black",
           position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = fill_colors, drop = FALSE) +
  scale_x_discrete(guide = offset_guide()) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.05)),
    guide = offset_guide()
  ) +
  labs(x = NULL, y = "Dependent Site Count", fill = "Tissue") +
  theme(
    axis.title.y = element_text(colour = "black", size = 22),
    # Keep X tick labels, but remove the X axis line and tick marks.
    axis.line.x = element_blank(),
    axis.ticks.x = element_blank(),
    # X tick labels are mutant gene names: render in Arial Bold Italic.
    axis.text.x = element_text(family = "Arial", colour = "black", size = 20,
                               angle = 45, hjust = 1, vjust = 1,
                               face = "bold.italic"),
    axis.text.y = element_text(colour = "black", size = 20),
    legend.text = element_text(family = "Arial", face = "bold", size = 18),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.key.height = unit(14, "points"),
    legend.margin = margin(6, 6, 6, 6)
  ) +
  coord_cartesian(clip = "off")

save_plot(p, file.path(script_dir, "dependent_site_count_stacked_barplot"), width = 10, height = 7)
message("[plot] Done: dependent_site_count_stacked_barplot.pdf/.png/.svg")
