#!/usr/bin/env Rscript
# chloroplast_chr_overview — 拟南芥叶绿体染色体（ChrC）环形 Ψ 位点总览
#
# 移植自 archive/20260902/paper/chloroplast_RNA/chloroplast_chr_overview/plot.R。
# 外圈 = Araport11_fixed.gtf 的 ChrC rRNA/tRNA/mRNA 基因；内圈径向 bar = WT leaf
# 位点的修饰比例。位点先按 Ch*/tRNA-Pt-*/ATCG* 映射到 ChrC genomic coordinate。
#
# 数据源：
#   rRNA  polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv   (chr ^Ch)
#   tRNA  tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv      (chr ^tRNA-Pt-)
#   mRNA  polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv   (chr ^ATCG)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(ggforce)
  library(ggprism)
  library(tibble)
})

feature_priority <- c(tRNA = 1L, rRNA = 2L, mRNA = 3L)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]
  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }
  normalizePath(getwd())
}

extract_attr <- function(attr, key) {
  pattern <- paste0(key, ' "([^"]+)"')
  stringr::str_match(attr, pattern)[, 2]
}

read_gtf <- function(path) {
  readr::read_tsv(
    path,
    comment = "#",
    col_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"),
    col_types = readr::cols(.default = readr::col_character(), start = readr::col_integer(), end = readr::col_integer()),
    show_col_types = FALSE
  ) %>%
    filter(seqname == "ChrC") %>%
    mutate(
      gene_id = extract_attr(attribute, "gene_id"),
      transcript_id = extract_attr(attribute, "transcript_id"),
      width = end - start + 1L
    )
}

read_site_input <- function(path, rna_type, chr_pattern) {
  dat <- readr::read_tsv(path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "ratio", "depth")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      site_strand = as.character(strand),
      ratio = as.numeric(ratio),
      depth = as.numeric(depth)
    ) %>%
    filter(grepl(chr_pattern, chr), is.finite(ratio), ratio >= 0, ratio <= 1) %>%
    transmute(
      rna_type = factor(rna_type, levels = c("rRNA", "tRNA", "mRNA")),
      chr,
      pos,
      site_strand,
      ratio,
      ratio_pct = ratio * 100,
      depth,
      site_id = paste(chr, pos, site_strand, sep = ":")
    )
}

run_analysis <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  register_font()

  theme_set(
    theme_base(base_size = SIZE_TOKENS$base_size, border = FALSE, sizes = SIZE_TOKENS) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA),
        panel.border = element_blank(),
        axis.title = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.margin = margin(0, 0, 0, 0)
      )
  )

  gtf_file <- file.path(repo_root, "Araport11_fixed.gtf")
  site_inputs <- list(
    rRNA = file.path(repo_root, "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv"),
    tRNA = file.path(repo_root, "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv"),
    mRNA = file.path(repo_root, "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv")
  )

  missing_inputs <- c(gtf_file, unlist(site_inputs))[!file.exists(c(gtf_file, unlist(site_inputs)))]
  if (length(missing_inputs) > 0) {
    stop("Missing required inputs:\n", paste(missing_inputs, collapse = "\n"))
  }

  gtf <- read_gtf(gtf_file)
  chr_len <- max(gtf$end, na.rm = TRUE)

  feature_tbl <- gtf %>%
    filter(feature %in% c("tRNA", "rRNA", "mRNA")) %>%
    transmute(
      chr = seqname,
      gene_id,
      transcript_id,
      feature_type = feature,
      start,
      end,
      strand,
      width,
      plot_label = gene_id,
      priority = feature_priority[feature]
    ) %>%
    arrange(start, end)

  if (nrow(feature_tbl) == 0) {
    stop("No ChrC tRNA/rRNA/mRNA features were found in ", gtf_file)
  }

  exon_tbl <- gtf %>%
    filter(feature == "exon", !is.na(transcript_id)) %>%
    select(transcript_id, exon_start_raw = start, exon_end_raw = end) %>%
    inner_join(
      feature_tbl %>%
        select(transcript_id, gene_id, feature_type, feature_start = start, feature_end = end, strand),
      by = "transcript_id"
    ) %>%
    transmute(
      transcript_id,
      gene_id,
      feature_type,
      feature_start,
      feature_end,
      strand,
      start = exon_start_raw,
      end = exon_end_raw
    ) %>%
    arrange(transcript_id, if_else(strand == "+", start, -end)) %>%
    group_by(transcript_id) %>%
    mutate(
      exon_len = end - start + 1L,
      transcript_offset = lag(cumsum(exon_len), default = 0L)
    ) %>%
    ungroup() %>%
    select(
      transcript_key = transcript_id, gene_id, feature_type,
      feature_start, feature_end, strand,
      exon_start = start, exon_end = end, exon_len, transcript_offset
    )

  gene_exon_tbl <- exon_tbl %>% mutate(transcript_key = gene_id)

  gene_map <- feature_tbl %>% select(gene_id, feature_type, start, end, strand)

  manual_alias_map <- tibble::tribble(
    ~chr, ~gene_id,
    "Ch16S", "ATCG00920",
    "Ch23S", "ATCG00950",
    "Ch4.5S", "ATCG00940",
    "Ch5S", "ATCG00930",
    "tRNA-Pt-Ala-TGC-1", "ATCG00200",
    "tRNA-Pt-Arg-ACG-1", "ATCG00710",
    "tRNA-Pt-Arg-TCT-1", "ATCG00110",
    "tRNA-Pt-Asn-GTT-1", "ATCG01140",
    "tRNA-Pt-Asp-GTC-1", "ATCG01100",
    "tRNA-Pt-Cys-GCA-1", "ATCG00840",
    "tRNA-Pt-eMet-CAT-1", "ATCG00970",
    "tRNA-Pt-Gln-TTG-1", "ATCG00620",
    "tRNA-Pt-Glu-TTC-1", "ATCG00370",
    "tRNA-Pt-Gly-GCC-1", "ATCG00260",
    "tRNA-Pt-Gly-TCC-1", "ATCG01150",
    "tRNA-Pt-His-GTG-1", "ATCG01160",
    "tRNA-Pt-Ile-CAT-1", "ATCG00610",
    "tRNA-Pt-Ile-GAT-1", "ATCG00090",
    "tRNA-Pt-iMet-CAT-1", "ATCG00820",
    "tRNA-Pt-Leu-CAA-1", "ATCG00630",
    "tRNA-Pt-Leu-TAA-1", "ATCG01040",
    "tRNA-Pt-Leu-TAG-1", "ATCG00410",
    "tRNA-Pt-Lys-TTT-1", "ATCG01030",
    "tRNA-Pt-Phe-GAA-1", "ATCG00010",
    "tRNA-Pt-Pro-TGG-1", "ATCG00460",
    "tRNA-Pt-Ser-CGA-1", "ATCG01250",
    "tRNA-Pt-Ser-GCT-1", "ATCG00030",
    "tRNA-Pt-Ser-GGA-1", "ATCG00700",
    "tRNA-Pt-Ser-TGA-1", "ATCG00100",
    "tRNA-Pt-Thr-GGT-1", "ATCG00830",
    "tRNA-Pt-Thr-TGT-1", "ATCG00880",
    "tRNA-Pt-Trp-CCA-1", "ATCG01020",
    "tRNA-Pt-Tyr-GTA-1", "ATCG01200",
    "tRNA-Pt-Val-GAC-1", "ATCG00590",
    "tRNA-Pt-Val-TAC-1", "ATCG00990"
  ) %>%
    left_join(gene_map, by = "gene_id") %>%
    select(chr, gene_id, feature_type, start, end, strand)

  site_tbl_raw <- bind_rows(
    read_site_input(site_inputs$rRNA, "rRNA", "^Ch"),
    read_site_input(site_inputs$tRNA, "tRNA", "^tRNA-Pt-"),
    read_site_input(site_inputs$mRNA, "mRNA", "^ATCG")
  ) %>%
    arrange(rna_type, chr, pos, site_strand)

  map_with_exons <- function(site_tbl, exon_map) {
    site_tbl %>%
      inner_join(exon_map, by = "transcript_key", relationship = "many-to-many") %>%
      filter(pos > transcript_offset, pos <= transcript_offset + exon_len) %>%
      mutate(
        exon_pos = pos - transcript_offset,
        chr_pos = case_when(
          strand == "+" ~ exon_start + exon_pos - 1L,
          strand == "-" ~ exon_end - exon_pos + 1L,
          TRUE ~ NA_integer_
        ),
        chr_pos = as.integer(chr_pos),
        source_chr = chr,
        mapped = !is.na(chr_pos) & chr_pos >= exon_start & chr_pos <= exon_end
      )
  }

  direct_sites <- site_tbl_raw %>%
    filter(str_detect(chr, "^ATCG")) %>%
    mutate(transcript_key = chr) %>%
    map_with_exons(exon_tbl)

  alias_sites <- site_tbl_raw %>%
    filter(!str_detect(chr, "^ATCG")) %>%
    left_join(manual_alias_map, by = "chr") %>%
    mutate(transcript_key = gene_id) %>%
    # 丢掉别名映射里的 gene_id：让 exon 表提供唯一一份，避免 join 后
    # gene_id.x/gene_id.y 把输出字段冲成 NA。
    select(-gene_id, -feature_type, -start, -end, -strand) %>%
    map_with_exons(gene_exon_tbl)

  site_tbl_mapped <- bind_rows(direct_sites, alias_sites) %>%
    distinct(rna_type, source_chr, gene_id, pos, chr_pos, site_strand, .keep_all = TRUE)

  unmapped_tbl <- site_tbl_raw %>%
    anti_join(
      site_tbl_mapped %>% transmute(rna_type, chr = source_chr, pos, site_strand),
      by = c("rna_type", "chr", "pos", "site_strand")
    )

  if (nrow(unmapped_tbl) > 0) {
    warning(nrow(unmapped_tbl), " chloroplast sites could not be mapped to ChrC coordinates and will be omitted.")
  }

  site_tbl <- site_tbl_mapped %>%
    filter(mapped) %>%
    select(rna_type, source_chr, gene_id, feature_type, pos, chr_pos, site_strand, ratio, ratio_pct, depth, site_id) %>%
    arrange(rna_type, chr_pos, source_chr, pos)

  if (nrow(site_tbl) == 0) {
    stop("No WT leaf chloroplast sites could be mapped to ChrC coordinates.")
  }

  feature_summary_tbl <- feature_tbl %>% count(feature_type, name = "n_features")
  site_summary_tbl <- site_tbl %>%
    group_by(rna_type) %>%
    summarise(n_sites = n(), mean_ratio_pct = mean(ratio_pct), .groups = "drop")

  readr::write_tsv(feature_tbl, file.path(script_dir, "chloroplast_chr_overview_features.tsv"))
  readr::write_tsv(site_tbl, file.path(script_dir, "chloroplast_chr_overview_sites.tsv"))
  readr::write_tsv(unmapped_tbl, file.path(script_dir, "chloroplast_chr_overview_unmapped_sites.tsv"))
  readr::write_tsv(feature_summary_tbl, file.path(script_dir, "chloroplast_chr_overview_feature_summary.tsv"))
  readr::write_tsv(site_summary_tbl, file.path(script_dir, "chloroplast_chr_overview_site_summary.tsv"))

  # ---------------------------------------------------------------------------
  # 环形布局
  # ---------------------------------------------------------------------------
  angle_offset <- -pi / 2
  pos_to_angle <- function(pos) angle_offset + 2 * pi * (pos - 1) / chr_len

  arc_points <- function(df, inner_radius, outer_radius, n = 72) {
    purrr::pmap_dfr(
      df,
      function(chr, gene_id, transcript_id, feature_type, start, end, strand, width, plot_label, priority, ...) {
        theta <- seq(pos_to_angle(start), pos_to_angle(end), length.out = n)
        tibble(
          gene_id = gene_id,
          transcript_id = transcript_id,
          feature_type = feature_type,
          start = start,
          end = end,
          strand = strand,
          order = seq_along(theta),
          x_inner = inner_radius * cos(theta),
          y_inner = inner_radius * sin(theta),
          x_outer = outer_radius * cos(theta),
          y_outer = outer_radius * sin(theta)
        )
      }
    )
  }

  feature_ring_inner_radius <- 0.62
  feature_ring_outer_radius <- 0.72

  feature_ring <- arc_points(feature_tbl, inner_radius = feature_ring_inner_radius, outer_radius = feature_ring_outer_radius) %>%
    group_by(gene_id, transcript_id, feature_type, start, end, strand) %>%
    reframe(
      x = c(x_inner, rev(x_outer)),
      y = c(y_inner, rev(y_outer)),
      vertex_order = seq_along(x)
    )

  chr_backbone <- tibble(theta = seq(0, 2 * pi, length.out = 720)) %>%
    mutate(x = feature_ring_outer_radius * cos(theta), y = feature_ring_outer_radius * sin(theta))

  site_point_radius <- 0.54
  site_bar_outer_radius <- 0.47
  site_bar_max_length <- 0.15

  site_bar_tbl <- site_tbl %>%
    mutate(
      theta = pos_to_angle(chr_pos),
      r0 = site_bar_outer_radius,
      r1 = r0 - site_bar_max_length * pmax(0, pmin(ratio, 1)),
      x = r0 * cos(theta),
      y = r0 * sin(theta),
      xend = r1 * cos(theta),
      yend = r1 * sin(theta)
    )

  site_point_tbl <- site_tbl %>%
    mutate(
      theta = pos_to_angle(chr_pos),
      r = site_point_radius,
      x = r * cos(theta),
      y = r * sin(theta)
    )

  mRNA_site_labels <- site_tbl %>%
    filter(rna_type == "mRNA", !is.na(gene_id)) %>%
    group_by(gene_id) %>%
    summarise(chr_pos = mean(chr_pos), .groups = "drop") %>%
    arrange(chr_pos, gene_id) %>%
    mutate(
      theta = pos_to_angle(chr_pos),
      label = gene_id,
      label_radius = rep(c(0.88, 1.05, 1.20), length.out = n()),
      x = label_radius * cos(theta),
      y = label_radius * sin(theta),
      xend = feature_ring_outer_radius * cos(theta),
      yend = feature_ring_outer_radius * sin(theta),
      angle_deg = theta * 180 / pi,
      text_angle = if_else(angle_deg > 90 & angle_deg < 270, angle_deg + 180, angle_deg),
      hjust = if_else(cos(theta) >= 0, 0, 1)
    )

  feature_colors <- c(rRNA = "#8DD3C7CC", mRNA = "#FDB462CC", tRNA = "#B3DE69CC")

  p <- ggplot() +
    geom_path(
      data = chr_backbone,
      aes(x = x, y = y),
      linewidth = 0.7,
      color = "black"
    ) +
    geom_polygon(
      data = feature_ring,
      aes(x = x, y = y, group = gene_id, fill = feature_type),
      color = "black",
      linewidth = 0.22
    ) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = feature_ring_inner_radius),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.55
    ) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = feature_ring_outer_radius),
      inherit.aes = FALSE,
      color = "black",
      linewidth = 0.55
    ) +
    geom_segment(
      data = site_bar_tbl,
      aes(x = x, y = y, xend = xend, yend = yend),
      color = "black",
      linewidth = 0.95,
      lineend = "round"
    ) +
    geom_segment(
      data = site_bar_tbl,
      aes(x = x, y = y, xend = xend, yend = yend, color = rna_type),
      linewidth = 0.55,
      lineend = "round"
    ) +
    geom_point(
      data = site_point_tbl,
      aes(x = x, y = y, fill = rna_type),
      shape = 21,
      color = "black",
      size = 1.75,
      stroke = 0.32,
      show.legend = FALSE
    ) +
    ggforce::geom_circle(
      aes(x0 = 0, y0 = 0, r = site_bar_outer_radius),
      color = "black",
      linewidth = 0.45,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = feature_colors, name = NULL, breaks = c("tRNA", "mRNA", "rRNA"), drop = FALSE) +
    scale_color_manual(values = feature_colors, guide = "none", drop = FALSE) +
    guides(fill = guide_legend(order = 1, override.aes = list(linewidth = 0.35))) +
    theme(
      legend.position = c(0.5, 0.5),
      legend.justification = c(0.5, 0.5),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 7),
      legend.key.size = unit(6, "points"),
      legend.spacing.y = unit(2, "points")
    )

  p_plain <- p +
    coord_equal(xlim = c(-0.78, 0.78), ylim = c(-0.78, 0.78), clip = "off")

  save_plot(p_plain, file.path(script_dir, "chloroplast_chr_overview"), width = 7.2, height = 7.2)

  p_annotated <- p +
    geom_segment(
      data = mRNA_site_labels,
      aes(x = xend, y = yend, xend = x, yend = y),
      color = "black",
      linewidth = 0.24,
      inherit.aes = FALSE
    ) +
    geom_text(
      data = mRNA_site_labels,
      aes(x = x, y = y, label = label, angle = text_angle, hjust = hjust),
      size = 2.3,
      color = "black",
      family = "Arial",
      inherit.aes = FALSE
    ) +
    coord_equal(xlim = c(-1.33, 1.33), ylim = c(-1.33, 1.33), clip = "off") +
    theme(plot.margin = margin(8, 8, 8, 8))

  save_plot(p_annotated, file.path(script_dir, "chloroplast_chr_overview_annotated"), width = 7.2, height = 7.2)
}

if (sys.nframe() == 0) {
  run_analysis()
}
