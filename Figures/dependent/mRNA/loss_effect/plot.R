#!/usr/bin/env Rscript
# ============================================================================
# plot.R — mRNA dependent 位点修饰下降效应 dotplot（7 突变体 × 3 组织，共 21 张）
#
# 输入：dependent_site_loss_effect_dotplot.tsv（data_prep.sh 生成）
#   每行一个 WT 阳性位点：WT_Ratio / Mut_Ratio / site_status(dependent|wt_only)
# 输出：<genotype>_<tissue>_dependent_site_loss_dotplot.pdf / .png / .svg
#   X 轴 = 突变体平均修饰比例，Y 轴 = WT 修饰比例；
#   灰点 = 全部 WT 阳性位点，红点 = 该基因型×组织的 dependent 位点。
#   坐标轴：分离式（offset axes，用户明确要求）。
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(ggtext)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

input_file <- "dependent_site_loss_effect_dotplot.tsv"
if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file, "\nRun data_prep.sh first.")
}

genotype_levels <- c("trua1a", "trua1b", "trua2", "trua3", "trua4", "trua5", "trua6")
tissue_levels <- c("Root", "leaf", "Flower")
ratio_breaks <- seq(0, 1, by = 0.2)
status_colors <- setNames(PAL_STATUS, c("WT sites", "Dependent sites"))

plot_df <- read.table(input_file, header = TRUE, sep = "\t", quote = "",
                      comment.char = "", stringsAsFactors = FALSE,
                      check.names = FALSE) |>
  mutate(
    genotype = factor(genotype, levels = genotype_levels),
    tissue = factor(tissue, levels = tissue_levels),
    WT_Ratio = as.numeric(WT_Ratio),
    Mut_Ratio = as.numeric(Mut_Ratio),
    Delta_Ratio = as.numeric(Delta_Ratio),
    site_status = factor(site_status,
                         levels = c("wt_only", "dependent"),
                         labels = c("WT sites", "Dependent sites"))
  )

if (nrow(plot_df) == 0) {
  stop("No rows found in plotting table: ", input_file)
}

# 用户明确要求分离坐标轴（offset axes）：
# theme_base_offset() 设置 axis.line + panel.border = blank
# 每个 scale 手动添加 guide_prism_offset()
init_theme_offset(base_size = 24)

# 单张面板：X = 突变体修饰比例，Y = WT 修饰比例，y=x 虚线参考
build_plot <- function(df, genotype) {
  df <- df |> arrange(site_status)   # dependent 排在最后绘制，压在灰点之上
  # ggtext markdown keeps the mutant gene italic while inheriting the global
  # Arial Bold style; the italic token is rendered as Arial Bold Italic.
  x_label <- sprintf("*%s* Modification Ratio", genotype)

  ggplot(df, aes(x = Mut_Ratio, y = WT_Ratio)) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.6,
                linetype = "dashed", color = "grey45") +
    geom_point(aes(fill = "WT sites"), size = 2.8, shape = 21,
               color = "black", stroke = 0.45) +
    geom_point(data = filter(df, site_status == "Dependent sites"),
               aes(fill = site_status), size = 3.0, shape = 21,
               color = "black", stroke = 0.55) +
    scale_fill_manual(values = status_colors,
                      breaks = c("WT sites", "Dependent sites"),
                      name = NULL) +
    scale_x_continuous(limits = c(0, 1), breaks = ratio_breaks,
                       labels = percent_format(accuracy = 1),
                       expand = expansion(mult = 0.02),
                       guide = ggprism::guide_prism_offset()) +
    scale_y_continuous(limits = c(0, 1), breaks = ratio_breaks,
                       labels = percent_format(accuracy = 1),
                       expand = expansion(mult = 0.02),
                       guide = ggprism::guide_prism_offset()) +
    coord_equal() +
    labs(x = x_label,
         y = "WT Modification Ratio") +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      legend.text = element_text(size = 14),
      legend.key = element_rect(fill = "white", colour = NA),
      axis.title.x = ggtext::element_markdown(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(10, 20, 10, 12)
    )
}

for (genotype in genotype_levels) {
  for (tissue in tissue_levels) {
    panel_df <- plot_df |>
      filter(.data$genotype == .env$genotype, .data$tissue == .env$tissue)
    if (nrow(panel_df) == 0) {
      stop("Missing plotting data for: ", genotype, " / ", tissue)
    }

    stem <- paste0(genotype, "_", tolower(tissue), "_dependent_site_loss_dotplot")
    save_plot(build_plot(panel_df, genotype), stem, width = 8, height = 6.04)
  }
}

message("[plot] Done: 21 dependent-site loss-effect dotplots with offset axes (pdf/png/svg)")
