#!/usr/bin/env Rscript
# ============================================================================
# plot.R — dependent 位点两方法 (methylKit vs dependent) overlap Venn 图
#   按组织分图：每个组织一张图，图内按 7 个基因型各画一个小 Venn。
#
# 输入：dependent_site_overlap.tsv  — tissue | genotype | method | site（data_prep.sh 生成）
# 输出：dependent_site_overlap_venn_<tissue>.pdf / .png / .svg
#       dependent_site_overlap_summary.tsv
# ============================================================================

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
})

# 关闭 VennDiagram 的 INFO 日志输出
if (requireNamespace("futile.logger", quietly = TRUE)) {
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

# ---- 唯一数据来源：data_prep.sh 生成的 TSV ----
df <- read.table("dependent_site_overlap.tsv", header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)

tissues <- c("Root", "leaf", "Flower")
tissue_display <- c(Root = "Roots", leaf = "Leaves", Flower = "Flowers")
genotypes <- c("trua1a", "trua1b", "trua2", "trua3", "trua4", "trua5", "trua6")

methods <- c("methylKit", "dependent")
fill <- PAL_METHOD

build_panel <- function(tissue, genotype) {
  sets <- lapply(methods, function(m) {
    unique(df$site[df$tissue == tissue & df$genotype == genotype & df$method == m])
  })
  names(sets) <- methods

  # 两方法均为 0：VennDiagram 不支持全空集合，用占位文本
  if (all(lengths(sets) == 0)) {
    return(grid::grobTree(
      grid::textGrob(genotype, x = 0.5, y = 0.82,
                     gp = grid::gpar(fontfamily = "Arial", fontface = "bold.italic", fontsize = 16)),
      grid::textGrob("0", x = 0.5, y = 0.40,
                     gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 20))
    ))
  }

  venn.diagram(
    x = sets,
    category.names = c("", ""),
    filename = NULL,
    disable.logging = TRUE,
    # 两个圆固定大小，不随集合大小/交集缩放
    scaled = FALSE,
    euler.d = FALSE,
    # 圆圈样式
    fill = as.character(fill),
    alpha = 0.45,
    col = "black",
    lwd = 1.2,
    lty = 1,
    # 数字
    cex = 1.7,
    fontface = "bold",
    fontfamily = "Arial",
    # 标题（基因型名，斜体）
    main = genotype,
    main.cex = 1.8,
    main.fontface = "bold.italic",
    main.fontfamily = "Arial",
    main.pos = c(0.5, 0.80),
    margin = 0.12
  )
}

# 底部图例（颜色注释）：红 = methylKit，蓝 = dependent
draw_legend <- function() {
  gp_text <- grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 16)
  sq <- 0.16  # 方块边长（英寸）
  grid::grid.rect(x = grid::unit(0.40, "npc"), y = grid::unit(0.05, "npc"),
                  width = grid::unit(sq, "inches"), height = grid::unit(sq, "inches"),
                  gp = grid::gpar(fill = PAL_METHOD[["methylKit"]], col = "black", lwd = 1))
  grid::grid.text("methylKit", x = grid::unit(0.418, "npc"), y = grid::unit(0.05, "npc"),
                  just = "left", gp = gp_text)
  grid::grid.rect(x = grid::unit(0.56, "npc"), y = grid::unit(0.05, "npc"),
                  width = grid::unit(sq, "inches"), height = grid::unit(sq, "inches"),
                  gp = grid::gpar(fill = PAL_METHOD[["dependent"]], col = "black", lwd = 1))
  grid::grid.text("dependent", x = grid::unit(0.578, "npc"), y = grid::unit(0.05, "npc"),
                  just = "left", gp = gp_text)
}

# 每个组织一张图：顶部组织标题 + 2x4 网格（7 个基因型，第 8 格留空）+ 底部图例
draw_tissue <- function(tissue) {
  panels <- lapply(genotypes, function(g) build_panel(tissue, g))

  grid::grid.newpage()
  grid::pushViewport(grid::viewport())
  grid::grid.text(tissue_display[[tissue]], x = 0.5, y = 0.97,
                  gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 22))
  # 正方形网格：respect=TRUE 使每个 cell 宽高比 1:1（圆保持正圆）
  grid::pushViewport(grid::viewport(
    x = grid::unit(0.5, "npc"),
    y = grid::unit(0.54, "npc"),
    width = grid::unit(11.2, "inches"),
    height = grid::unit(5.6, "inches"),
    layout = grid::grid.layout(2, 4, respect = TRUE)
  ))
  for (i in seq_along(panels)) {
    row <- (i - 1) %/% 4 + 1
    col <- (i - 1) %% 4 + 1
    grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
    grid::grid.draw(panels[[i]])
    grid::popViewport()
  }
  grid::popViewport()
  # 底部图例
  draw_legend()
  grid::popViewport()
}

# ---- 保存 PDF / PNG / SVG（每个组织一张）----
for (tissue in tissues) {
  save_plot(
    plot = function() draw_tissue(tissue),  # draw_tissue 自带 grid.newpage()
    filename = paste0("dependent_site_overlap_venn_", tissue),
    width = 12,
    height = 7
  )
}

# ---- 汇总 TSV ----
summaries <- lapply(tissues, function(tissue) {
  do.call(rbind, lapply(genotypes, function(genotype) {
    s <- lapply(methods, function(m) {
      unique(df$site[df$tissue == tissue & df$genotype == genotype & df$method == m])
    })
    names(s) <- methods
    n <- lengths(s)
    inter <- length(intersect(s[["methylKit"]], s[["dependent"]]))
    data.frame(
      tissue          = tissue,
      genotype        = genotype,
      methylKit_total = n[["methylKit"]],
      dependent_total = n[["dependent"]],
      overlap         = inter,
      methylKit_only  = n[["methylKit"]] - inter,
      dependent_only  = n[["dependent"]] - inter
    )
  }))
})
summary_df <- do.call(rbind, summaries)
write.table(summary_df, "dependent_site_overlap_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("[venn] Done: dependent_site_overlap_venn_{Root,leaf,Flower}.pdf/.png/.svg")
print(summary_df, row.names = FALSE)
