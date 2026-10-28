#!/usr/bin/env Rscript
# plot.R — Sprinzl 39 位点图组（Flower 组；一张热图 + 一张修饰比例分布图）
#
# 输入（跑 data_prep.sh 生成）：pos39_flower_site_table.tsv
# 实现全在 ../site_specificity_common.R，本文件只声明参数。
# 图面组织按项目统一口径显示 Flowers；底层目录/数据仍用 Flower。

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
common_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
shared_dir <- normalizePath(file.path(common_dir, "..", "..", "..", "shared"), mustWork = TRUE)

source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))
source(file.path(common_dir, "site_specificity_common.R"))

genotypes <- c("trua3", "trua6")
sprinzl_pos <- 39L
position_label <- paste0("Ψ", sprinzl_pos)
input_path <- file.path(script_dir, "pos39_flower_site_table.tsv")
file_prefix <- "pos39_flower"

# ---- 热图：行 = WT Ψ39 修饰的核代表 tRNA（WT 比例降序），列 = WT + 各株系 ----
heatmap_result <- run_site_specificity_heatmap(
  input_path = input_path,
  genotypes = genotypes,
  out_dir = script_dir,
  file_prefix = file_prefix
)
print(heatmap_result$metadata)

# ---- 修饰比例分布：x = WT / trua3_1 / trua3_2 / trua6_1 / trua6_2 ----
distribution_result <- run_site_specificity_distribution(
  input_path = input_path,
  genotypes = genotypes,
  position_label = position_label,
  out_dir = script_dir,
  file_prefix = file_prefix
)
print(distribution_result$summary)
