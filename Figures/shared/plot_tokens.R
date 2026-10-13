# plot_tokens.R — Figures 共享绘图常量(尺寸与调色板)
# 纯数据文件:零依赖、零副作用,直接 source 即可。
# 机制(字体/主题/导出)在 plot_style.R;两文件互不依赖。
#
# 约定:颜色用 #RRGGBBCC(alpha = 0.8),除非个别历史值本身不带 alpha。

# ============================================================================
# 1. 尺寸与线宽(wildtype_landscape 显式字号体系)
# ============================================================================
SIZE_TOKENS <- list(
  base_linewidth = 0.6,
  base_size = 24,
  axis_text_size = 20,
  axis_title_size = 22,
  annotation_text_size = 16,
  legend_text_size = 18,
  legend_title_size = 18,
  strip_text_size = 18,
  plot_title_size = 18
)

# ============================================================================
# 2. 调色板(仅收录跨脚本复用 >= 2 处的色板;单用色板留在各 plot.R 本地)
# ============================================================================

# 组织三色,粉彩版(dependent/dependent_methylKit 的 sites_count、ratio_distribution)
# 不命名:调用方名字不统一(Root/Leaf/Flower 与 Roots/Leaves/Flowers 两种),
# 各脚本自行 setNames()
PAL_TISSUE <- c("#FF9FB2CC", "#7ED6B5CC", "#8FB8FFCC")

# 组织三色,NPG 不透明版(保留备用)
PAL_TISSUE_NPG <- c(
  "Root" = "#E64B35FF",
  "Leaf" = "#00A087FF",
  "Flower" = "#3C5488FF"
)

# 位点状态二色:非 dependent/WT(灰)、dependent(红)
# 不命名:调用方叫法不一("WT sites" / "non-dependent sites"),各脚本自行 setNames()
PAL_STATUS <- c("#B8B8B8CC", "#E64B35CC")

# mRNA 区域三色(dependent 两树的 region_distribution)
PAL_REGION <- c(
  "5'UTR" = "#F6D55CCC",
  "CDS" = "#F4A7B9CC",
  "3'UTR" = "#8FD3F4CC"
)

# 方法对比二色(dependent_methylKit 的 sites_count/overlap Venn;原值不带 alpha)
PAL_METHOD <- c(
  "methylKit" = "#E64B35",
  "dependent" = "#3C5488"
)
