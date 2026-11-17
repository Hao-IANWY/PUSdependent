# sci-plot theme.R — 出版级 R/ggplot2 默认主题
# 所有生成的 R 脚本通过 source("theme.R") 引用此文件，保证风格一致
#
# 铁律：R 脚本只通过 read.table() 读取 TSV 数据，不内嵌数据。
#       本文件不读取任何数据，只定义主题、字体、颜色和输出设备。

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(sysfonts)
  library(svglite)
})

# ============================================================================
# 1. Arial 字体配置
# ============================================================================
ARIAL_FONT_PATH <- "/home/HL/miniforge3/envs/analysis/fonts/arial.ttf"

if (!file.exists(ARIAL_FONT_PATH)) {
  stop(
    "Default Arial font file not found at ",
    ARIAL_FONT_PATH, "\n",
    "Please locate an Arial font file with: fc-list | grep -i 'arial'"
  )
}

if (!"Arial" %in% sysfonts::font_families()) {
  sysfonts::font_add("Arial", regular = ARIAL_FONT_PATH)
}

# ============================================================================
# 2. PDF 输出设备（可编辑 Arial 文本，禁止 showtext）
# ============================================================================
cairo_pdf_arial <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename = filename,
    width    = width,
    height   = height,
    family   = "Arial",
    ...
  )
}

# ============================================================================
# 2b. SVG 输出设备（使用 svglite，文本保持为完整 <text> 元素）
# ============================================================================
svg_arial <- function(filename, width, height, ...) {
  svglite::svglite(
    filename = filename,
    width    = width,
    height   = height,
    fix_text_size = FALSE,
    ...
  )
}

# ============================================================================
# 3. 默认主题 — non-negotiable rules
# ============================================================================
# base_size 的 16 只是兜底缺省值，不要直接吃。字号由画布宽度决定：
#     BASE_SIZE <- round(FIG_WIDTH_IN * 72 / 27)   # 即 画布宽(in) × 2.67
#     9 in -> 24pt | 6 in -> 16pt | 12 in -> 32pt
# 因为图总是按固定宽度呈现（报告 / 栏宽 / 幻灯片），只有「字 : 画布宽」这个
# 比值决定字看起来多大，画布绝对英寸数不影响。详见 SKILL.md
# 「字号必须由画布宽度决定」。
theme_sciplot <- function(base_size = 16, border = TRUE) {
  ggprism::theme_prism(base_size = base_size, border = border) +
    theme(
      # 全局：Arial、bold、黑色
      text             = element_text(family = "Arial", face = "bold", colour = "black"),
      # 无网格
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 黑色面板边框，linewidth 0.6
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      # 图例：可见边框，与面板边框一致
      legend.background     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
      # 离散图例间距
      legend.key.spacing.y  = unit(12, "points"),
      # 坐标轴刻度：与面板边框一致
      axis.ticks        = element_line(color = "black", linewidth = 0.6)
    )
}

# 应用默认主题
theme_set(theme_sciplot())

# ============================================================================
# 4. Offset Axes 变体（偏移/分离坐标轴）
# ============================================================================
theme_sciplot_offset <- function(base_size = 16) {
  ggprism::theme_prism(base_size = base_size, border = FALSE) +
    theme(
      text             = element_text(family = "Arial", face = "bold", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      legend.background     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.key.spacing.y  = unit(12, "points"),
      axis.ticks        = element_line(color = "black", linewidth = 0.6)
    )
}

# 应用 offset axes 主题的便捷函数
use_offset_theme <- function(base_size = 16) {
  theme_set(theme_sciplot_offset(base_size))
}

# 为 scale_* 添加 offset guide 的辅助函数
offset_guide <- function() {
  ggprism::guide_prism_offset()
}

# ============================================================================
# 5. Y 轴刻度上限 — 最高刻度必须罩住最高柱（不许柱子捅出刻度）
# ============================================================================
# ggplot 默认由「带留白的面板范围」反推刻度，于是同一个 expand 在不同量级的
# 数据上会翻出两种车：柱子捅出最高刻度（刻度 75 而柱子 94），或最高刻度远高于
# 柱子（刻度 300 而柱子 253，面板一大片空白）。
# 这里把刻度改成显式给定：pretty() 会把范围外扩到整数刻度，其最大值必然 >=
# 数据最大值；再把面板上界正好收到该刻度处，刻度与留白解耦，两种翻车同时消失。
#
# 用法：
#   ax <- y_axis_cover(bar_totals[["total"]])   # 传所有会被画出来的 y 值
#   scale_y_continuous(breaks = ax$breaks, expand = ax$expand)
#
# values: 堆叠柱状图必须传各柱的累计总量（柱顶），不是单段值，否则面板会算矮。
# n:      目标刻度个数，pretty() 会自行取整（默认 4）。
y_axis_cover <- function(values, n = 4, lower_mult = 0.02) {
  values <- values[is.finite(values)]
  lo <- min(c(0, values))
  hi <- max(c(0, values))

  # 退化情形（全 0 / 空数据）：给一根 0-1 的安全轴
  if (hi <= lo) {
    return(list(breaks = c(0, 1), top = 1,
                expand = ggplot2::expansion(mult = c(lower_mult, lower_mult))))
  }

  breaks <- pretty(c(lo, hi), n = n)
  top    <- max(breaks)

  list(
    breaks = breaks,
    top    = top,
    # 上侧留白由“最高刻度 - 最高柱”决定，面板正好收到最高刻度
    expand = ggplot2::expansion(mult = c(lower_mult, (top - hi) / (hi - lo)))
  )
}

# ============================================================================
# 6. 默认调色板
# ============================================================================

# 糖果/粉彩调色板（candy/pastel，首选用色）
# 使用 #RRGGBBCC 格式：alpha = 0.8 (CC)
PAL_CANDY <- c(
  "#E64B35CC",  # 红
  "#00A087CC",  # 绿
  "#3C5488CC",  # 蓝
  "#F39B7FCC",  # 橙
  "#8491B4CC",  # 紫灰
  "#91D1C2CC",  # 青
  "#DC0000CC",  # 深红
  "#7E6148CC"   # 棕
)

# 组织/组织类型专用色
COLOR_ROOT   <- "#E64B35CC"   # 红
COLOR_LEAF   <- "#00A087CC"   # 绿
COLOR_FLOWER <- "#3C5488CC"   # 蓝

# Volcano plot 专用色
COLOR_UP        <- "#E64B35CC"   # 上调：红
COLOR_DOWN      <- "#3C5488CC"   # 下调：蓝
COLOR_NS        <- "#B0B0B0CC"   # 不显著：灰
COLOR_SPIKE     <- "#00A087CC"   # spike-in：绿

# Okabe-Ito 色盲安全 8 色
PAL_OKABE <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# 热度图：白→红渐变
HEATMAP_GRADIENT <- c("white", "#E64B35CC")

# ============================================================================
# 7. 通用保存函数 — 同时输出 PDF、PNG 和 SVG
# ============================================================================
save_sciplot <- function(plot, filename, width = 8, height = 5,
                         formats = c("pdf", "png", "svg")) {
  base_name <- tools::file_path_sans_ext(filename)
  ext <- tools::file_ext(filename)
  if (ext != "") {
    base_name <- tools::file_path_sans_ext(filename)
  } else {
    base_name <- filename
  }

  if ("pdf" %in% formats) {
    pdf_file <- paste0(base_name, ".pdf")
    ggsave(
      filename = pdf_file,
      plot     = plot,
      width    = width,
      height   = height,
      device   = cairo_pdf_arial
    )
    message("[sci-plot] Wrote: ", pdf_file)
  }
  if ("png" %in% formats) {
    png_file <- paste0(base_name, ".png")
    ggsave(
      filename = png_file,
      plot     = plot,
      width    = width,
      height   = height,
      dpi      = 600,
      type     = "cairo"
    )
    message("[sci-plot] Wrote: ", png_file)
  }
  if ("svg" %in% formats) {
    svg_file <- paste0(base_name, ".svg")
    svg_arial(svg_file, width = width, height = height)
    print(plot)
    grDevices::dev.off()
    message("[sci-plot] Wrote: ", svg_file)
  }
}

save_sciplot_all <- save_sciplot  # 向后兼容别名

message("[sci-plot] theme.R loaded — Arial, ggprism, cairo_pdf ready")
