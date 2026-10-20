#!/usr/bin/env Rscript
# plot.R — TRUA1A-dependent rRNA motif sequence logos with offset axes and Bits scores
# 输入：trua1a_motif_sequences.tsv（data_prep.sh 生成）
# 输出：有位点组织的 trua1a.<tissue>.dependent.seqlogo.pdf/.png/.svg
# 与 mRNA 版的差异：rRNA 位点集中在 Flower（Root/leaf 常为 0），空组织跳过而非报错。

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggseqlogo)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))

input_file <- file.path(script_dir, "trua1a_motif_sequences.tsv")
if (!file.exists(input_file)) {
  stop("Input table not found: ", input_file, "\nRun data_prep.sh first.")
}

expected_tissues <- c(Root = "Roots", leaf = "Leaves", Flower = "Flowers")
export_width <- 6
export_height <- 5

motif_df <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c("genotype", "tissue", "chr", "pos", "strand", "sequence")
missing_columns <- setdiff(required_columns, names(motif_df))
if (length(missing_columns) > 0) {
  stop("Required columns missing from ", input_file, ": ",
       paste(missing_columns, collapse = ", "))
}
if (nrow(motif_df) == 0) {
  stop("No sequences found in ", input_file)
}
if (any(motif_df$genotype != "trua1a")) {
  stop("Unexpected genotype in ", input_file, "; this workflow is scoped to trua1a.")
}
if (any(!motif_df$tissue %in% names(expected_tissues))) {
  stop("Unexpected tissue in ", input_file, ": ",
       paste(sort(unique(motif_df$tissue[!motif_df$tissue %in% names(expected_tissues)])),
             collapse = ", "))
}
if (any(nchar(motif_df$sequence) != 5L) || any(grepl("[^ACGT]", motif_df$sequence))) {
  stop("All sequences must be 5-nt DNA contexts containing only A, C, G, and T.")
}
if (anyDuplicated(motif_df[, c("tissue", "chr", "pos", "strand")]) != 0) {
  stop("Duplicate tissue/site keys found in ", input_file)
}

tissues_with_sites <- names(expected_tissues)[
  vapply(names(expected_tissues), function(t) sum(motif_df$tissue == t) > 0L, logical(1))
]
if (length(tissues_with_sites) == 0L) {
  stop("No sequences found for any tissue in ", input_file)
}
skipped_tissues <- setdiff(names(expected_tissues), tissues_with_sites)
if (length(skipped_tissues) > 0) {
  message("[plot] Skipping tissues without dependent sites: ",
          paste(skipped_tissues, collapse = ", "))
}

plot_logo <- function(seqs, tissue_label, output_stem) {
  p <- ggseqlogo::ggseqlogo(
    seqs,
    method = "bits",
    seq_type = "dna",
    col_scheme = "nucleotide"
  ) +
    scale_x_continuous(
      breaks = 1:5,
      labels = c("-2", "-1", "0", "1", "2"),
      expand = expansion(mult = c(0.08, 0.08)),
      guide = ggprism::guide_prism_offset()
    ) +
    scale_y_continuous(
      limits = c(0, 2),
      breaks = seq(0, 2, by = 0.5),
      expand = expansion(mult = c(0.04, 0.08)),
      guide = ggprism::guide_prism_offset()
    ) +
    labs(title = tissue_label, x = NULL, y = "Bits") +
    theme(
      legend.position = "none",
      axis.title.y = element_text(size = 22, face = "bold", margin = margin(r = 10)),
      axis.title.x = element_blank(),
      axis.text.x = element_text(face = "bold", size = 20),
      axis.text.y = element_text(face = "bold", size = 20),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.ticks.x = element_line(color = "black", linewidth = 0.6),
      axis.ticks.y = element_line(color = "black", linewidth = 0.6),
      axis.ticks.length = unit(3, "mm"),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      plot.background = element_rect(color = NA, fill = "white", linewidth = 0),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18,
                                margin = margin(b = 8)),
      plot.margin = margin(20, 20, 18, 18)
    ) +
    coord_cartesian(clip = "off")

  # save_plot() 只剥离已知图形扩展名(pdf/png/svg),".seqlogo" 词干可直传
  save_plot(
    p,
    file.path(script_dir, output_stem),
    width = export_width,
    height = export_height
  )
}

init_theme_offset(base_size = 24)

for (tissue in tissues_with_sites) {
  seqs <- motif_df$sequence[motif_df$tissue == tissue]
  output_stem <- paste0("trua1a.", tissue, ".dependent.seqlogo")
  plot_logo(seqs, unname(expected_tissues[[tissue]]), output_stem)
}

message("[plot] Done: TRUA1A-dependent rRNA motif sequence logos (pdf/png/svg)")
