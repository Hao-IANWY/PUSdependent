#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(trackViewer)
  library(GenomicRanges)
  library(grid)
})

plot_width <- 16
plot_height <- 6
y_scale_factor <- 0.05

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]

  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(active_path)) {
      return(dirname(normalizePath(active_path)))
    }
  }

  normalizePath(getwd())
}

setup_style <- function(repo_root) {
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  invisible("Arial")
}

read_positive_sites <- function(file_path, tissue_name) {
  readr::read_tsv(file_path, show_col_types = FALSE) |>
    transmute(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      Tissue = tissue_name,
      ratio = as.numeric(ratio),
      depth = as.numeric(depth)
    )
}

make_transcript_track <- function(chr, chr_length, fill, feature_height = unit(4.5, "mm")) {
  features <- GRanges(
    seqnames = chr,
    ranges = IRanges(1, chr_length, names = chr)
  )
  features$fill <- fill
  features$border <- "#8A8A8A"
  features$height <- feature_height
  features$featureLayerID <- chr
  features
}

make_site_track <- function(site_df, chr) {
  chr_df <- site_df |>
    filter(chr == !!chr) |>
    arrange(pos)

  fill_colors <- dplyr::case_when(
    chr_df$ratio * 100 < 20 ~ "#8FD3E8",
    chr_df$ratio * 100 <= 50 ~ "#F2B36E",
    TRUE ~ "#F28B82"
  )

  gr <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = chr_df$pos, width = 1, names = paste0(chr, "_", chr_df$pos)),
    score = pmax(chr_df$ratio * 100 * y_scale_factor, 1)
  )
  gr$color <- fill_colors
  gr$border <- "#4D4D4D"
  gr$alpha <- 0.95
  gr$cex <- 1.0
  gr$node.label <- sprintf("%d", round(chr_df$ratio * 100, 0))
  gr$node.label.col <- "#222222"
  gr$node.label.cex <- 0.75
  gr$node.label.fontface <- "bold"
  gr$label.parameter.rot <- 90
  gr$label.parameter.col <- fill_colors
  gr$label.parameter.label <- names(gr)
  gr$label.parameter.cex <- 0.72
  gr$label.parameter.draw <- TRUE
  gr
}

build_xaxis <- function(chr_length) {
  xaxis_at <- unique(c(1, pretty(c(1, chr_length), n = 4), chr_length))
  xaxis_at <- xaxis_at[xaxis_at >= 1 & xaxis_at <= chr_length]
  xaxis_at <- sort(as.integer(round(xaxis_at)))
  names(xaxis_at) <- xaxis_at
  xaxis_at
}

draw_lolliplot <- function(site_track, transcript_track, chr_length, show_yaxis = TRUE, show_ylab = TRUE) {
  xaxis_at <- build_xaxis(chr_length)
  yaxis_at <- c(0, 20, 40, 60, 80, 100)
  yaxis_pos <- yaxis_at * y_scale_factor
  names(yaxis_pos) <- yaxis_at

  trackViewer::lolliplot(
    site_track,
    features = transcript_track,
    ranges = GRanges(as.character(seqnames(transcript_track))[1], IRanges(1, chr_length)),
    type = "circle",
    newpage = FALSE,
    ylab = if (show_ylab) expression(Psi ~ "ratio (%)") else "",
    ylab.gp = gpar(col = "black", fontsize = SIZE_TOKENS$axis_title_size, fontface = "bold", fontfamily = "Arial"),
    yaxis = if (show_yaxis) yaxis_pos else FALSE,
    yaxis.gp = gpar(col = "black", fontsize = SIZE_TOKENS$axis_text_size, fontface = "bold", fontfamily = "Arial"),
    xaxis = xaxis_at,
    xaxis.gp = gpar(col = "#5A5A5A", fontsize = SIZE_TOKENS$axis_text_size, fontface = "bold", fontfamily = "Arial"),
    legend = NULL,
    cex = 1.15,
    rescale = FALSE,
    label_on_feature = FALSE,
    lollipop_style_switch_limit = 1000
  )
}

draw_transcript_panel <- function(site_df, chr, chr_length, transcript_fill, show_yaxis = TRUE, show_ylab = TRUE) {
  transcript_df <- site_df |>
    filter(chr == !!chr) |>
    arrange(pos)

  if (nrow(transcript_df) == 0) {
    grid::grid.text(
      label = paste(chr, "(no sites)"),
      x = 0.5,
      y = 0.5,
      gp = gpar(col = "#4D4D4D", fontsize = SIZE_TOKENS$annotation_text_size, fontface = "bold", fontfamily = "Arial")
    )
    return(invisible(NULL))
  }

  site_track <- make_site_track(transcript_df, chr = chr)
  transcript_track <- make_transcript_track(chr = chr, chr_length = chr_length, fill = transcript_fill)

  draw_lolliplot(
    site_track = site_track,
    transcript_track = transcript_track,
    chr_length = chr_length,
    show_yaxis = show_yaxis,
    show_ylab = show_ylab
  )
}

draw_rrna_layout <- function(site_df, transcript_lengths, transcript_fill_colors) {
  panel_order <- c("25S", "18S", "5.8S")
  panel_widths <- sqrt(as.numeric(transcript_lengths[panel_order]))

  grid::pushViewport(
    grid::viewport(
      layout = grid::grid.layout(
        nrow = 3,
        ncol = 3,
        heights = unit.c(
          unit(0.9, "null"),
          unit(1, "mm"),
          unit(0.9, "null")
        ),
        widths = unit.c(
          unit(panel_widths[2], "null"),
          unit(6, "mm"),
          unit(panel_widths[3], "null")
        )
      )
    )
  )

  grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1:3))
  draw_transcript_panel(
    site_df = site_df,
    chr = "25S",
    chr_length = transcript_lengths[["25S"]],
    transcript_fill = transcript_fill_colors[["25S"]],
    show_yaxis = TRUE,
    show_ylab = TRUE
  )
  grid::upViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
  draw_transcript_panel(
    site_df = site_df,
    chr = "18S",
    chr_length = transcript_lengths[["18S"]],
    transcript_fill = transcript_fill_colors[["18S"]],
    show_yaxis = TRUE,
    show_ylab = FALSE
  )
  grid::upViewport()

  grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 3))
  draw_transcript_panel(
    site_df = site_df,
    chr = "5.8S",
    chr_length = transcript_lengths[["5.8S"]],
    transcript_fill = transcript_fill_colors[["5.8S"]],
    show_yaxis = FALSE,
    show_ylab = FALSE
  )
  grid::upViewport()

  grid::upViewport()
}

run_single_tissue_rrna_distribution <- function(
  tissue = "Flower",
  out_dir = NULL,
  output_stem = NULL,
  save_outputs = TRUE
) {
  tissue_options <- c("Root", "Leaf", "Flower")
  tissue <- as.character(tissue)[1]
  if (!tissue %in% tissue_options) {
    stop("Unsupported tissue: ", tissue, ". Supported tissues: ", paste(tissue_options, collapse = ", "))
  }

  script_dir <- get_script_dir()
  project_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = TRUE)
  setup_style(project_root)
  out_dir <- if (is.null(out_dir)) script_dir else normalizePath(out_dir, mustWork = FALSE)

  input_path <- switch(
    tissue,
    Root = file.path(project_root, "polyARNA/STAT/fisher/rRNA/WT-Root_sites.positive.tsv"),
    Leaf = file.path(project_root, "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv"),
    Flower = file.path(project_root, "polyARNA/STAT/fisher/rRNA/WT-Flower_sites.positive.tsv")
  )

  transcript_lengths <- c("18S" = 1808, "5.8S" = 164, "25S" = 3375)
  transcript_fill_colors <- c(
    "25S" = "#FCE8E6",
    "18S" = "#E6F0FA",
    "5.8S" = "#FDF1DC"
  )

  if (is.null(output_stem)) {
    output_stem <- paste0("single_tissue_rrna_distribution_", tolower(tissue))
  }

  site_df <- read_positive_sites(input_path, tissue) |>
    filter(chr %in% names(transcript_lengths))

  if (nrow(site_df) == 0) {
    stop("No positive rRNA sites found for tissue: ", tissue)
  }
  outputs <- list(
    plot_pdf = file.path(out_dir, paste0(output_stem, ".pdf")),
    plot_svg = file.path(out_dir, paste0(output_stem, ".svg")),
    plot_png = file.path(out_dir, paste0(output_stem, ".png"))
  )

  if (save_outputs) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    save_plot(
      plot = function() {
        grid::grid.newpage()
        draw_rrna_layout(
          site_df = site_df,
          transcript_lengths = transcript_lengths,
          transcript_fill_colors = transcript_fill_colors
        )
      },
      filename = file.path(out_dir, output_stem),
      width = plot_width,
      height = plot_height
    )
  }

  list(
    tissue = tissue,
    input_path = input_path,
    plot = NULL,
    panel_grob = NULL,
    output_paths = outputs
  )
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  tissue <- if (length(args) >= 1) args[1] else "Flower"
  invisible(run_single_tissue_rrna_distribution(tissue = tissue))
}
