#!/usr/bin/env Rscript
# chloroplast_site_overview — WT leaf 叶绿体 rRNA / tRNA / mRNA Ψ 位点总览
#
# 移植自 archive/20260902/paper/chloroplast_RNA/chloroplast_site_overview/plot.R。
# 数据源换成当前流程：
#   rRNA  polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv   (chr == Ch23S)
#   tRNA  tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv      (chr ^tRNA-Pt-)
#   mRNA  polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv   (chr ^ATCG)
# 每个 RNA 类型单独一张修饰比例密度图。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

region_fill <- c(
  "rRNA" = "#8DD3C7CC",
  "tRNA" = "#B3DE69CC",
  "mRNA" = "#FDB462CC"
)

rna_titles <- c(
  rRNA = "Chloroplast rRNA",
  tRNA = "Chloroplast tRNA",
  mRNA = "Chloroplast mRNA"
)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]
  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }
  normalizePath(getwd())
}

read_chloroplast_sites <- function(file_path, rna_type, chr_pattern) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "ratio", "depth")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      ratio = as.numeric(ratio),
      depth = as.numeric(depth)
    ) %>%
    filter(grepl(chr_pattern, chr)) %>%
    transmute(
      rna_type = factor(rna_type, levels = c("rRNA", "tRNA", "mRNA")),
      chr,
      pos,
      strand,
      ratio,
      depth,
      ratio_pct = ratio * 100,
      site_id = paste(chr, pos, strand, sep = ":")
    ) %>%
    filter(is.finite(ratio_pct), ratio_pct >= 0, ratio_pct <= 100)
}

compute_peak_density <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2 || length(unique(x)) < 2) {
    return(0.02)
  }
  max(stats::density(x, from = 0, to = 100, n = 512)$y, na.rm = TRUE)
}

run_analysis <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)

  input_specs <- list(
    list(
      file = file.path(repo_root, "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv"),
      rna_type = "rRNA",
      chr_pattern = "^Ch23S$"
    ),
    list(
      file = file.path(repo_root, "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv"),
      rna_type = "tRNA",
      chr_pattern = "^tRNA-Pt-"
    ),
    list(
      file = file.path(repo_root, "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv"),
      rna_type = "mRNA",
      chr_pattern = "^ATCG"
    )
  )

  missing_inputs <- vapply(input_specs, function(x) !file.exists(x$file), logical(1))
  if (any(missing_inputs)) {
    stop(
      "Missing required inputs:\n",
      paste(vapply(input_specs[missing_inputs], `[[`, character(1), "file"), collapse = "\n")
    )
  }

  site_tbl <- bind_rows(lapply(input_specs, function(spec) {
    read_chloroplast_sites(spec$file, spec$rna_type, spec$chr_pattern)
  })) %>%
    arrange(rna_type, chr, pos, strand)

  if (nrow(site_tbl) == 0) {
    stop("No chloroplast sites were found for the requested WT leaf filters.")
  }

  count_tbl <- site_tbl %>%
    count(rna_type, name = "n_sites") %>%
    mutate(label = paste0("n = ", scales::comma(n_sites)))

  readr::write_tsv(site_tbl, file.path(script_dir, "chloroplast_site_overview_density_sites.tsv"))
  readr::write_tsv(count_tbl, file.path(script_dir, "chloroplast_site_overview_count_summary.tsv"))

  label_tbl <- site_tbl %>%
    group_by(rna_type) %>%
    summarise(peak_density = compute_peak_density(ratio_pct), .groups = "drop") %>%
    left_join(count_tbl, by = "rna_type") %>%
    mutate(x = 50, y = peak_density * 0.95)

  build_rna_plot <- function(rna_name) {
    rna_df <- site_tbl %>% filter(rna_type == rna_name)
    rna_label <- label_tbl %>% filter(rna_type == rna_name)
    ax <- y_axis_cover(c(0, rna_label$peak_density))

    ggplot(rna_df, aes(x = ratio_pct, fill = rna_type)) +
      geom_density(
        color = "black",
        linewidth = SIZE_TOKENS$base_linewidth,
        alpha = 0.7,
        trim = TRUE
      ) +
      geom_text(
        data = rna_label,
        aes(x = x, y = y, label = label),
        inherit.aes = FALSE,
        family = "Arial",
        fontface = "bold",
        size = SIZE_TOKENS$annotation_text_size / ggplot2::.pt,
        hjust = 0.5,
        vjust = 1
      ) +
      scale_fill_manual(values = region_fill, drop = FALSE) +
      scale_x_continuous(
        limits = c(0, 100),
        breaks = seq(0, 100, by = 20),
        expand = expansion(mult = c(0, 0)),
        name = "Modification ratio (%)"
      ) +
      scale_y_continuous(
        name = "Density",
        breaks = ax$breaks,
        expand = ax$expand
      ) +
      labs(title = rna_titles[[as.character(rna_name)]]) +
      coord_cartesian(clip = "off") +
      theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = SIZE_TOKENS$plot_title_size)
      )
  }

  for (rna_name in c("rRNA", "tRNA", "mRNA")) {
    stem <- file.path(
      script_dir,
      paste0("chloroplast_site_overview_density_", tolower(rna_name))
    )
    save_plot(build_rna_plot(rna_name), stem, width = 8, height = 6.6)
  }
}

if (sys.nframe() == 0) {
  run_analysis()
}
