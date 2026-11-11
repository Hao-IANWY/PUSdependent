#!/usr/bin/env Rscript
# Sequence logo of WT mRNA Ψ-site windows (y axis = bits; each column of a
# 5-mer peaks at 2 bits).
# Panels: Roots / Leaves / Flowers / All (three tissues pooled), styled to
# match the dependent motif pipeline of the repository
# (Figures/Tissue/dependent/motif/plot_dependent_motif_seqlogo_panel.R):
# Arial bold, theme_prism, offset guides, and the four-color ACGT scheme.
# Unlike the reference script, which uses method = "prob" (percent y axis),
# this figure uses the bits scale.
#
# Usage: Rscript plot_mrna_motif_seqlogo.R <sequences.tsv> <output_dir> <window_width>
#   window_width = total window length (±2 -> 5, ±5 -> 11, ±10 -> 21)

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(ggseqlogo)
  library(scales)
  library(sysfonts)
  library(patchwork)
  library(svglite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript plot_mrna_motif_seqlogo.R <sequences.tsv> <output_dir> <window_width>")
}
input_tsv <- args[[1]]
output_dir <- args[[2]]
window_width <- as.integer(args[[3]])
if (is.na(window_width) || window_width < 5 || window_width %% 2 != 1) {
  stop("window_width must be an odd integer >= 5, got: ", args[[3]])
}

if (!file.exists(input_tsv)) stop("Input file not found: ", input_tsv)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Arial font (consistent with the other figures in this repository)
arial_font_path <- "/home/HL/miniforge3/envs/analysis/fonts/arial.ttf"
if (!file.exists(arial_font_path)) {
  stop("Arial font not found at ", arial_font_path)
}
if (!"Arial" %in% sysfonts::font_families()) {
  sysfonts::font_add("Arial", regular = arial_font_path)
}

cairo_pdf_arial <- function(filename, width, height, ...) {
  grDevices::cairo_pdf(filename = filename, width = width, height = height, family = "Arial", ...)
}
svg_arial <- function(filename, width, height, ...) {
  svglite::svglite(filename = filename, width = width, height = height, fix_text_size = FALSE, ...)
}

df <- read.delim(input_tsv, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
required_cols <- c("Tissue", "sequence")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) stop("Input missing columns: ", paste(missing_cols, collapse = ", "))

tissue_levels <- c("Root", "Leaf", "Flower", "All")
df$Tissue <- factor(df$Tissue, levels = tissue_levels)
df$sequence <- toupper(gsub("U", "T", df$sequence))
width_pattern <- paste0("^[ACGT]{", window_width, "}$")
df <- df[!is.na(df$Tissue) & grepl(width_pattern, df$sequence), ]
if (nrow(df) == 0) stop("No usable motif sequences found.")

# Pooled panel: all sequences from the three tissues
all_df <- rbind(
  df[df$Tissue %in% c("Root", "Leaf", "Flower"), ],
  transform(df[df$Tissue %in% c("Root", "Leaf", "Flower"), ], Tissue = factor("All", levels = tissue_levels))
)
df <- rbind(df[df$Tissue %in% c("Root", "Leaf", "Flower"), ], all_df[all_df$Tissue == "All", ])

sequence_widths <- unique(nchar(df$sequence))
if (length(sequence_widths) != 1 || sequence_widths != window_width) {
  stop("Sequence widths inconsistent with window_width: ",
       paste(sequence_widths, collapse = ","))
}
half_width <- window_width %/% 2
position_labels <- seq_len(window_width) - (half_width + 1)  # -half .. +half
# x-axis tick throttling: label every position for ±2, every 2 nt for ±5,
# and every 5 nt for ±10 and wider, to avoid crowding
x_labels <- if (half_width <= 2) {
  position_labels
} else if (half_width <= 5) {
  unique(seq(-half_width, half_width, 2))
} else {
  unique(c(seq(-half_width, half_width, 5), 0))
}
x_breaks <- match(x_labels, position_labels)

# Four-color ACGT scheme (matching the dependent-motif pipeline)
base_colors <- make_col_scheme(
  chars = c("A", "C", "G", "T"),
  groups = c("A", "C", "G", "T"),
  cols = c("#109648CC", "#255C99CC", "#F7B32BCC", "#D62839CC")
)

panel_labels <- c(Root = "Roots", Leaf = "Leaves", Flower = "Flowers", All = "All tissues")

plot_list <- list()
for (tissue in tissue_levels) {
  seqs <- df$sequence[df$Tissue == tissue]
  is_first <- (tissue == tissue_levels[1])

  p <- ggseqlogo(seqs, method = "bits", seq_type = "dna", col_scheme = base_colors) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(add = 0.22),
      guide = ggprism::guide_prism_offset()
    ) +
    scale_y_continuous(
      limits = c(0, 2.05),
      breaks = seq(0, 2, by = 0.5),
      expand = expansion(mult = c(0.02, 0.04)),
      guide = ggprism::guide_prism_offset()
    ) +
    labs(
      x = NULL,
      y = "Bits",
      title = panel_labels[[tissue]],
      subtitle = sprintf("(n=%d)", length(seqs))
    ) +
    theme_prism(base_size = 22, border = FALSE) +
    theme(
      text = element_text(family = "Arial", face = "bold", colour = "black"),
      plot.title = element_text(family = "Arial", face = "bold.italic", hjust = 0.5, size = 22,
                                margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 16, margin = margin(b = 6)),
      plot.margin = margin(12, 6, 4, 6),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black", linewidth = 1.0),
      axis.ticks = element_line(color = "black", linewidth = 1.0),
      axis.ticks.length = unit(5, "pt"),
      legend.position = "none"
    ) +
    coord_cartesian(clip = "off")

  plot_list[[tissue]] <- p
  cat(sprintf("Created panel: %s (n=%d)\n", tissue, length(seqs)))
}

combined_plot <- wrap_plots(plot_list, nrow = 1) +
  patchwork::plot_annotation(
    caption = "Position relative to Ψ site (nt)",
    theme = theme(
      plot.caption = element_text(
        hjust = 0.5, family = "Arial", face = "bold", colour = "black", size = 20
      ),
      plot.margin = margin(4, 8, 8, 8)
    )
  )

output_width <- 15
output_height <- 5.2

cat("Saving combined seqlogo figure...\n")

ggsave(
  file.path(output_dir, "mrna_motif_seqlogo.pdf"),
  combined_plot, width = output_width, height = output_height, device = cairo_pdf_arial
)
ggsave(
  file.path(output_dir, "mrna_motif_seqlogo.png"),
  combined_plot, width = output_width, height = output_height, dpi = 600
)
svg_arial(file.path(output_dir, "mrna_motif_seqlogo.svg"), width = output_width, height = output_height)
print(combined_plot)
grDevices::dev.off()

cat("Figures saved:\n")
cat(sprintf("  - %s/mrna_motif_seqlogo.pdf\n", output_dir))
cat(sprintf("  - %s/mrna_motif_seqlogo.png\n", output_dir))
cat(sprintf("  - %s/mrna_motif_seqlogo.svg\n", output_dir))
