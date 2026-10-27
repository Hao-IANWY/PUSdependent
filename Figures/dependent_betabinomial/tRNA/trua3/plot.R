#!/usr/bin/env Rscript
# plot.R — 本基因型的 tRNA TRUA-dependent 位点图组
#
# 输入（跑 data_prep.sh 生成）：
#   <genotype>_all_tissues.dependent.tsv  组织拆分的 dependent 位点（含 fix_pos 与每株系修饰比）
#   <genotype>_all_tissues.full.dependent.tsv  .full. 版：Root 行替换为另两组织并集、修饰比取 Root 实测
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

genotype <- "trua3"

# ---- 组织拆分热图：行 = 组织 × 代表 tRNA，值 = Delta Ratio(WT - Mut)。7 个基因型都出 ----
tissue_result <- run_tissue_split_heatmap(
  genotype = genotype,
  base_dir = tRNA_dir,
  out_dir = script_dir
)
print(tissue_result$metadata)

source(file.path(tRNA_dir, "dependent_trna_distribution_common.R"))

# ---- Ψ 位点分布箱线图：x = 组织，y = 修饰比例(%)，WT / <genotype>_1 / <genotype>_2 三组 ----
distribution_result <- run_single_panel_dependent_position_distributions(
  genotype = genotype,
  input_path = file.path(script_dir, paste0(genotype, "_all_tissues.dependent.tsv")),
  out_dir = script_dir,
  plot_specs = tibble::tibble(
    compartment = c("Cy"),
    label = c(paste0("cy-tRNA \u03a8", "38/39/40")),
    fix_pos_values = list(c(38L, 39L, 40L)),
    output_suffix = c("cy_38_39_40")
  )
)
print(distribution_result$metadata)

# ---- .full. 版：Roots 组展示另外两个组织 dependent 位点的并集（修饰比取 Root 实测）；
#      Leaves / Flowers 组与常规版相同。数据由 data_prep.sh 生成 *all_tissues.full* 表 ----
full_distribution_result <- run_single_panel_dependent_position_distributions(
  genotype = genotype,
  input_path = file.path(script_dir, paste0(genotype, "_all_tissues.full.dependent.tsv")),
  out_dir = script_dir,
  file_tag = ".full",
  plot_specs = tibble::tibble(
    compartment = c("Cy"),
    label = c(paste0("cy-tRNA Ψ", "38/39/40")),
    fix_pos_values = list(c(38L, 39L, 40L)),
    output_suffix = c("cy_38_39_40.full")
  )
)
print(full_distribution_result$metadata)

source(file.path(tRNA_dir, "dependent_trna_motif_common.R"))

# ---- 位点 motif seqlogo：union 去重位点按 compartment × Sprinzl 位号选 5-nt 上下文。
#      motif 只取 38/39（老板指定）；trua3 数据里本就没有 Cy_40 位点，与 38/39/40 等价 ----
motif_result <- run_single_panel_dependent_motifs(
  genotype = genotype,
  input_path = file.path(script_dir, paste0(genotype, "_union.dependent.tsv")),
  out_dir = script_dir,
  plot_specs = tibble::tibble(
    compartment = c("Cy"),
    label = c(paste0("cy-tRNA Ψ", "38/39")),
    fix_pos_values = list(c(38L, 39L)),
    output_suffix = c("cy_38_39")
  )
)
print(motif_result$metadata)
