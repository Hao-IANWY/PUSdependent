# plot_style.R — Figures 共享绘图机制(字体 / 主题 / 导出)
# 纯机制文件:零依赖、零顶层副作用(source 不改变任何全局状态)。
# 尺寸与调色板常量在 plot_tokens.R;两文件互不依赖,可单独 source。
#
# 字体策略:不使用 showtext;PDF(cairo_pdf)、SVG(svglite)、
# PNG(png(type = "cairo")) 全部经 fontconfig 解析 "Arial"。
# register_font() 只做字体文件存在性检查 + sysfonts 注册(备 showtext 场景)。

# ============================================================================
# 1. Arial 字体
# ============================================================================
ARIAL_FONT_PATH <- "/home/HL/miniforge3/envs/analysis/fonts/arial.ttf"

register_font <- function() {
  if (!file.exists(ARIAL_FONT_PATH)) {
    stop(
      "Default Arial font file not found at ", ARIAL_FONT_PATH, "\n",
      "Please locate an Arial font file with: fc-list | grep -i 'arial'"
    )
  }
  suppressPackageStartupMessages(library(sysfonts))
  if (!"Arial" %in% sysfonts::font_families()) {
    sysfonts::font_add("Arial", regular = ARIAL_FONT_PATH)
  }
  invisible("Arial")
}

# ============================================================================
# 2. 输出设备
# ============================================================================
# PDF:可编辑 Arial 文本
cairo_pdf_arial <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename = filename,
    width    = width,
    height   = height,
    family   = "Arial",
    ...
  )
}

# SVG:文本保持为完整 <text> 元素
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
# 3. 主题
# ============================================================================
# sizes = NULL:纯 ggprism 比例字号(默认,原 sci-plot 行为)
# sizes = SIZE_TOKENS:显式字号体系(原 wildtype_landscape fig1 行为),
#                     SIZE_TOKENS 由 plot_tokens.R 提供
theme_base <- function(base_size = 16, border = TRUE, sizes = NULL) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggprism)
  })

  if (is.null(sizes)) {
    return(
      ggprism::theme_prism(base_size = base_size, border = border) +
        ggplot2::theme(
          # 全局:Arial、bold、黑色
          text             = element_text(family = "Arial", face = "bold", colour = "black"),
          # 无网格
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # 黑色面板边框,linewidth 0.6
          panel.border     = if (border) {
            element_rect(color = "black", fill = NA, linewidth = 0.6)
          } else {
            element_blank()
          },
          # 图例:可见边框,与面板边框一致
          legend.background     = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
          # 离散图例间距
          legend.key.spacing.y  = grid::unit(12, "points"),
          # 坐标轴刻度:与面板边框一致
          axis.ticks        = element_line(color = "black", linewidth = 0.6)
        )
    )
  }

  ggprism::theme_prism(base_size = base_size, border = border) +
    ggplot2::theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = sizes$axis_text_size, colour = "black"),
      axis.title = element_text(size = sizes$axis_title_size, colour = "black"),
      panel.border = if (border) {
        element_rect(color = "black", fill = NA, linewidth = sizes$base_linewidth)
      } else {
        element_blank()
      },
      axis.ticks = element_line(color = "black", linewidth = sizes$base_linewidth),
      legend.background = element_rect(color = "black", fill = NA, linewidth = sizes$base_linewidth),
      legend.box.background = element_rect(color = "black", fill = NA, linewidth = sizes$base_linewidth),
      legend.key.spacing.y = grid::unit(12, "points"),
      legend.text = element_text(size = sizes$legend_text_size, colour = "black"),
      legend.title = element_text(size = sizes$legend_title_size, colour = "black"),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      strip.background = element_rect(fill = "white", colour = "black", linewidth = sizes$base_linewidth),
      strip.text = element_text(size = sizes$strip_text_size, colour = "black"),
      plot.title = element_text(size = sizes$plot_title_size, hjust = 0.5, margin = margin(b = 6)),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8))
    )
}

# Offset Axes 变体(偏移/分离坐标轴)
theme_base_offset <- function(base_size = 16) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggprism)
  })

  ggprism::theme_prism(base_size = base_size, border = FALSE) +
    ggplot2::theme(
      text             = element_text(family = "Arial", face = "bold", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border     = element_blank(),
      axis.line        = element_line(color = "black", linewidth = 0.6),
      legend.background     = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.key.spacing.y  = grid::unit(12, "points"),
      axis.ticks        = element_line(color = "black", linewidth = 0.6)
    )
}

# 应用主题(显式 theme_set + 字体注册;本文件 source 时无任何副作用)
init_theme <- function(base_size = 16, border = TRUE, sizes = NULL) {
  suppressPackageStartupMessages(library(ggplot2))
  theme_set(theme_base(base_size = base_size, border = border, sizes = sizes))
  register_font()
  invisible("Arial")
}

init_theme_offset <- function(base_size = 16) {
  suppressPackageStartupMessages(library(ggplot2))
  theme_set(theme_base_offset(base_size = base_size))
  register_font()
  invisible("Arial")
}

# 为 scale_* 添加 offset guide 的辅助函数
offset_guide <- function() {
  ggprism::guide_prism_offset()
}

# ============================================================================
# 4. Y 轴刻度上限 — 最高刻度必须罩住最高的柱子
# ============================================================================
# ggplot 默认由「带留白的面板范围」反推刻度,于是同一个 expand 会在不同量级的
# 数据上翻出两种车:柱子捅出最高刻度(刻度 75 而柱子 94),或最高刻度远高于柱子
# (刻度 300 而柱子 253,面板一大片空白)。
# 这里把刻度改成显式给定:pretty() 会把范围外扩到整数刻度,其最大值必然 >=
# 数据最大值;再把面板上界正好收到该刻度处,刻度与留白解耦,两种翻车同时消失。
#
# 用法:
#   ax <- y_axis_cover(bar_totals[["total"]])   # 传所有会被画出来的 y 值
#   scale_y_continuous(breaks = ax$breaks, expand = ax$expand)
#
# values: 堆叠柱状图必须传各柱的累计总量(柱顶),不是单段值,否则面板会算矮。
# n:      目标刻度个数,pretty() 会自行取整(默认 4)。
y_axis_cover <- function(values, n = 4, lower_mult = 0.02) {
  values <- values[is.finite(values)]
  lo <- min(c(0, values))
  hi <- max(c(0, values))

  # 退化情形(全 0 / 空数据):给一根 0-1 的安全轴
  if (hi <= lo) {
    return(list(breaks = c(0, 1), top = 1,
                expand = ggplot2::expansion(mult = c(lower_mult, lower_mult))))
  }

  breaks <- pretty(c(lo, hi), n = n)
  top    <- max(breaks)

  list(
    breaks = breaks,
    top    = top,
    # 上侧留白由"最高刻度 - 最高柱"决定,面板正好收到最高刻度
    expand = ggplot2::expansion(mult = c(lower_mult, (top - hi) / (hi - lo)))
  )
}

# ============================================================================
# 5. 通用保存函数 — 同时输出 PDF、PNG 和 SVG
# ============================================================================
# plot 支持三种输入:
#   - ggplot 对象:走 ggsave
#   - grid grob:grid.newpage() + grid.draw()
#   - 零参绘制函数:直接调用(函数须自己在当前设备完成全部绘制,
#     如需清页自行 grid.newpage();适用于 viewport 过程式绘制)
# filename 语义:结尾是已知图形扩展名(pdf/png/svg)则剥掉,否则整体视为词干
#   ("x.pdf" 与 "x" 等价;"trua1a.Root.dependent.seqlogo" 原样保留)
save_plot <- function(plot, filename, width = 8, height = 5,
                      formats = c("pdf", "png", "svg"), dpi = 600) {
  suppressPackageStartupMessages(library(ggplot2))

  ext <- tolower(tools::file_ext(filename))
  stem <- if (ext %in% c("pdf", "png", "svg")) tools::file_path_sans_ext(filename) else filename

  draw <- if (is.function(plot)) {
    plot
  } else if (grid::is.grob(plot)) {
    function() {
      grid::grid.newpage()
      grid::grid.draw(plot)
    }
  } else {
    NULL
  }

  if ("pdf" %in% formats) {
    pdf_file <- paste0(stem, ".pdf")
    if (is.null(draw)) {
      ggsave(filename = pdf_file, plot = plot, width = width, height = height,
             device = cairo_pdf_arial)
    } else {
      cairo_pdf_arial(pdf_file, width = width, height = height)
      draw()
      grDevices::dev.off()
    }
    message("[plot_style] Wrote: ", pdf_file)
  }
  if ("png" %in% formats) {
    png_file <- paste0(stem, ".png")
    if (is.null(draw)) {
      ggsave(filename = png_file, plot = plot, width = width, height = height,
             dpi = dpi, type = "cairo")
    } else {
      # family = "Arial":兜住 grid 绘制中 fontfamily = "" 的文字(如 VennDiagram
      # 标签、trackViewer 节点数字),否则 cairo 默认 "sans" 会渲染成 DejaVu,
      # 与 PDF(cairo_pdf family = "Arial")不一致
      grDevices::png(png_file, width = width, height = height,
                     units = "in", res = dpi, type = "cairo", bg = "white",
                     family = "Arial")
      draw()
      grDevices::dev.off()
    }
    message("[plot_style] Wrote: ", png_file)
  }
  if ("svg" %in% formats) {
    svg_file <- paste0(stem, ".svg")
    svg_arial(svg_file, width = width, height = height)
    if (is.null(draw)) print(plot) else draw()
    grDevices::dev.off()
    message("[plot_style] Wrote: ", svg_file)
  }

  invisible(stem)
}
