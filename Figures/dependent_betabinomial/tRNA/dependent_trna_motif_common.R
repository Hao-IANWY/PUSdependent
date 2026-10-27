#!/usr/bin/env Rscript
# dependent_trna_motif_common.R — 单个基因型的 dependent 位点 motif seqlogo
#
# 与 dependent_trna_distribution_common.R 同一套选点口径（compartment × fix_pos 组合，
# plot_specs 三列同构），只是位点集换读 <genotype>_union.dependent.tsv（跨组织去重，
# 同一位点不因多组织重复而抬高 bits），上下文取 5-nt context（中心 0 位 = 修饰位点）。
#
# 每个 panel 一张 seqlogo：x = 相对修饰位点的偏移 -2..+2，y = bits（0–2），
# 样式与 ../motif/plot.R（基因型×组织全景版）一致。

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggprism)
  library(ggseqlogo)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

normalize_motif_specs <- function(plot_specs) {
  plot_specs |>
    mutate(
      compartment = as.character(compartment),
      label = as.character(label),
      fix_pos_values = purrr::map(fix_pos_values, as.integer)
    )
}

prepare_dependent_motif_data <- function(site_df, plot_specs) {
  required_cols <- c("chr", "pos", "fix_pos", "strand", "context")
  missing_cols <- setdiff(required_cols, names(site_df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in dependent-site table: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_specs <- normalize_motif_specs(plot_specs)

  base_df <- site_df |>
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      fix_pos = as.integer(fix_pos),
      compartment = dplyr::case_when(
        stringr::str_starts(chr, "tRNA-Nu-") ~ "Cy",
        stringr::str_starts(chr, "tRNA-Mt-") ~ "Mt",
        stringr::str_starts(chr, "tRNA-Pt-") ~ "Pt",
        TRUE ~ NA_character_
      ),
      # context 清洗与 motif/data_prep.sh 同口径：只留 ACGTU、5 长度、U 归一为 T
      sequence = toupper(as.character(context)),
      sequence = stringr::str_replace_all(sequence, "[^ACGTU]", "")
    ) |>
    filter(!is.na(compartment), stringr::str_length(sequence) == 5L) |>
    # 防御性去重：union 表本身每位点一行，换喂 all_tissues 时也不让跨组织重复计数
    dplyr::distinct(chr, pos, strand, .keep_all = TRUE) |>
    mutate(sequence = stringr::str_replace_all(sequence, "U", "T"))

  filtered_df <- purrr::map_dfr(seq_len(nrow(plot_specs)), function(i) {
    spec <- plot_specs[i, ]
    group_df <- base_df |>
      filter(
        compartment == spec$compartment[[1]],
        fix_pos %in% spec$fix_pos_values[[1]]
      ) |>
      mutate(position_group = spec$label[[1]])
    if (nrow(group_df) == 0) {
      stop(
        "No dependent sites matched compartment=", spec$compartment[[1]],
        " fix_pos in {", paste(spec$fix_pos_values[[1]], collapse = ","),
        "} — check the spec for ", spec$label[[1]], "."
      )
    }
    group_df
  })

  filtered_df |>
    dplyr::select(position_group, compartment, fix_pos, chr, pos, strand, sequence) |>
    arrange(position_group, chr, pos)
}

build_motif_logo_plot <- function(seqs, title_label) {
  ggseqlogo::ggseqlogo(
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
    labs(title = title_label, x = NULL, y = "Bits") +
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
}

run_single_panel_dependent_motifs <- function(
  genotype,
  input_path,
  out_dir,
  plot_specs,
  plot_width = 6,
  plot_height = 5
) {
  if (!file.exists(input_path)) {
    stop("Missing input table: ", input_path, "\n先跑 data_prep.sh 生成。")
  }

  plot_specs <- normalize_motif_specs(plot_specs)
  plot_df <- prepare_dependent_motif_data(
    site_df = readr::read_tsv(input_path, show_col_types = FALSE),
    plot_specs = plot_specs
  )

  init_theme_offset(base_size = 24)

  counts_df <- plot_df |>
    group_by(position_group) |>
    summarise(
      compartment = dplyr::first(compartment),
      fix_positions = paste(sort(unique(fix_pos)), collapse = "/"),
      n_sites = dplyr::n(),
      .groups = "drop"
    )

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  counts_tsv <- file.path(out_dir, paste0(genotype, "_dependent_motif_site_counts.tsv"))
  readr::write_tsv(counts_df, counts_tsv)

  output_tbl <- purrr::pmap_dfr(
    list(plot_specs$label, plot_specs$output_suffix),
    function(position_label, output_suffix) {
      seqs <- plot_df$sequence[plot_df$position_group == position_label]
      output_stem <- file.path(
        out_dir, paste0(genotype, "_", output_suffix, "_dependent_motif")
      )
      save_plot(
        plot = build_motif_logo_plot(seqs, position_label),
        filename = output_stem,
        width = plot_width,
        height = plot_height,
        formats = c("pdf", "png", "svg"),
        dpi = export_dpi
      )
      tibble::tibble(
        position_group = position_label,
        n_sites = length(seqs),
        plot_pdf = paste0(output_stem, ".pdf"),
        plot_png = paste0(output_stem, ".png"),
        plot_svg = paste0(output_stem, ".svg")
      )
    }
  )

  list(
    plot_data = plot_df,
    counts = counts_df,
    metadata = output_tbl,
    counts_tsv = counts_tsv
  )
}
