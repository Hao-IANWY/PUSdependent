# plot_style.R — shared plotting machinery for Figures (fonts / theme / export)
# Pure machinery file: zero dependencies and zero top-level side effects
# (sourcing it changes no global state).
# Size and palette constants live in plot_tokens.R; the two files are
# independent of each other and can be sourced separately.
#
# Font strategy: showtext is not used; PDF (cairo_pdf), SVG (svglite), and
# PNG (png(type = "cairo")) all resolve "Arial" through fontconfig.
# register_font() only verifies that the font file exists and registers it
# with sysfonts (kept for showtext-based scenarios).

# ============================================================================
# 1. Arial font
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
# 2. Output devices
# ============================================================================
# PDF: keeps text editable in Arial
cairo_pdf_arial <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(
    filename = filename,
    width    = width,
    height   = height,
    family   = "Arial",
    ...
  )
}

# SVG: text is preserved as intact <text> elements
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
# 3. Theme
# ============================================================================
# sizes = NULL: pure ggprism-scaled font sizes (default, classic sci-plot behavior)
# sizes = SIZE_TOKENS: explicit font-size system (the wildtype_landscape fig-1
#                      convention); SIZE_TOKENS is provided by plot_tokens.R
theme_base <- function(base_size = 16, border = TRUE, sizes = NULL) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ggprism)
  })

  if (is.null(sizes)) {
    return(
      ggprism::theme_prism(base_size = base_size, border = border) +
        ggplot2::theme(
          # Global: Arial, bold, black
          text             = element_text(family = "Arial", face = "bold", colour = "black"),
          # No grid lines
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          # Black panel border, linewidth 0.6
          panel.border     = if (border) {
            element_rect(color = "black", fill = NA, linewidth = 0.6)
          } else {
            element_blank()
          },
          # Legend: visible border matching the panel border
          legend.background     = element_rect(color = "black", fill = NA, linewidth = 0.6),
          legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
          # Key spacing for discrete legends
          legend.key.spacing.y  = grid::unit(12, "points"),
          # Axis ticks: match the panel border
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

# Offset-axes variant (shifted / detached axes)
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

# Apply the theme (explicit theme_set + font registration; sourcing this file
# has no side effects)
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

# Helper that attaches an offset guide to a scale_*
offset_guide <- function() {
  ggprism::guide_prism_offset()
}

# ============================================================================
# 4. Y-axis top tick — the top tick must cover the tallest bar
# ============================================================================
# By default ggplot derives breaks from the expanded panel range, so one fixed
# expand setting fails in two different ways depending on data magnitude:
# bars poke past the top tick (tick 75 while the bar reaches 94), or the top
# tick sits far above the bars (tick 300 for a bar of 253, leaving a large
# blank panel). Here breaks are supplied explicitly instead: pretty() rounds
# the range outward, so its maximum is necessarily >= the data maximum; the
# panel upper limit is then clipped exactly to that tick, decoupling ticks
# from the expansion and eliminating both failure modes.
#
# Usage:
#   ax <- y_axis_cover(bar_totals[["total"]])   # pass every y value that will be drawn
#   scale_y_continuous(breaks = ax$breaks, expand = ax$expand)
#
# values: for stacked bars pass each bar's cumulative total (the bar top),
#         not per-segment values, or the panel will come out too short.
# n:      target number of ticks; pretty() rounds as needed (default 4).
y_axis_cover <- function(values, n = 4, lower_mult = 0.02) {
  values <- values[is.finite(values)]
  lo <- min(c(0, values))
  hi <- max(c(0, values))

  # Degenerate case (all zero / empty data): return a safe 0-1 axis
  if (hi <= lo) {
    return(list(breaks = c(0, 1), top = 1,
                expand = ggplot2::expansion(mult = c(lower_mult, lower_mult))))
  }

  breaks <- pretty(c(lo, hi), n = n)
  # pretty() only guarantees "nice" breaks, not coverage of the data: on
  # small-magnitude values max(breaks) < hi can occur (e.g. a 0.0208 peak
  # yields breaks only up to 0.020). Append one extra tick, otherwise the
  # negative expansion would clip the top of the data.
  if (max(breaks) < hi) {
    breaks <- c(breaks, max(breaks) + diff(breaks)[length(breaks) - 1L])
  }
  top    <- max(breaks)

  # Pad the panel upper limit by "top tick - tallest bar" so it lands exactly
  # on the top tick. Add a tiny epsilon: floating-point error can leave the
  # panel limit 1e-18 below `top`, and ggplot then drops the tick sitting on
  # the boundary (the top tick disappears).
  upper_mult <- (top - hi) / (hi - lo)
  if (!is.finite(upper_mult) || upper_mult <= 0) upper_mult <- 0.05

  list(
    breaks = breaks,
    top    = top,
    expand = ggplot2::expansion(mult = c(lower_mult, upper_mult * (1 + 1e-4) + 1e-9))
  )
}

# ============================================================================
# 5. Generic save helper — write PDF, PNG, and SVG in one call
# ============================================================================
# `plot` accepts three kinds of input:
#   - a ggplot object: routed through ggsave
#   - a grid grob: grid.newpage() + grid.draw()
#   - a zero-argument drawing function: called directly (the function must
#     complete all drawing on the current device itself, calling
#     grid.newpage() first if a fresh page is needed; intended for procedural
#     viewport drawing)
# filename semantics: a trailing known graphics extension (pdf/png/svg) is
#   stripped; otherwise the whole string is used as the stem
#   ("x.pdf" is equivalent to "x"; "trua1a.Root.dependent.seqlogo" is kept as is)
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
      # family = "Arial": catch text drawn by grid with fontfamily = ""
      # (e.g. VennDiagram labels, trackViewer node counts); otherwise cairo's
      # default "sans" renders DejaVu Sans, inconsistent with the PDF
      # (cairo_pdf with family = "Arial")
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
