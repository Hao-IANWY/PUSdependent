#!/usr/bin/env Rscript
# plot.R — 本基因型的 tRNA TRUA-dependent 位点图组
#
# 输入（跑 data_prep.sh 生成）：
#   <genotype>_all_tissues.dependent.tsv  组织拆分的 dependent 位点（含 fix_pos 与每株系修饰比）
#   <genotype>_union.dependent.tsv        跨组织并集（仅 union 热图用）
# 实现全在 ../dependent_trna_heatmap_common.R 与 ../dependent_trna_distribution_common.R，
# 本文件只声明基因型、勾选本基因型要出的图。

suppressPackageStartupMessages({
  library(tibble)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grepl("^--file=", args)
  if (any(file_arg)) return(dirname(normalizePath(sub("^--file=", "", args[file_arg][1]))))
  normalizePath(getwd())
}

script_dir <- get_script_dir()
tRNA_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
shared_dir <- normalizePath(file.path(tRNA_dir, "..", "..", "shared"), mustWork = TRUE)

source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))
source(file.path(tRNA_dir, "dependent_trna_heatmap_common.R"))

genotype <- "trua5"

# ---- 组织拆分热图：行 = 组织 × 代表 tRNA，值 = Delta Ratio(WT - Mut)。7 个基因型都出 ----
tissue_result <- run_tissue_split_heatmap(
  genotype = genotype,
  base_dir = tRNA_dir,
  out_dir = script_dir
)
print(tissue_result$metadata)

# ---- 跨组织并集热图：行 = 代表 tRNA，值 = 该位点跨组织最大的 Delta Ratio ----
union_result <- run_union_dependent_heatmap(
  genotype = genotype,
  base_dir = tRNA_dir,
  out_dir = script_dir
)
print(union_result$metadata)

source(file.path(tRNA_dir, "dependent_trna_distribution_common.R"))

# ---- Ψ 位点分布箱线图：x = 组织，y = 修饰比例(%)，WT / <genotype>_1 / <genotype>_2 三组 ----
# 只取 canonical 28（anticodon stem 5' 端第一位）。trua5 在这一位没有 Ψ27 的
# dependent 位点，故不写成 trua1a 那样的 Ψ27/28 合并面板。
distribution_result <- run_single_panel_dependent_position_distributions(
  genotype = genotype,
  input_path = file.path(script_dir, paste0(genotype, "_all_tissues.dependent.tsv")),
  out_dir = script_dir,
  plot_specs = tibble::tibble(
    compartment = c("Cy"),
    label = c(paste0("cy-tRNA Ψ", "28")),
    fix_pos_values = list(28L),
    output_suffix = c("cy_28")
  )
)
print(distribution_result$metadata)
