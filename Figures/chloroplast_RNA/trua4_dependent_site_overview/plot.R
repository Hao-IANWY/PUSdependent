#!/usr/bin/env Rscript
# trua4_dependent_site_overview — 叶绿体全部 WT Ψ 位点在 WT vs trua4 中的修饰比例
#
# 移植自 archive/20260902/paper/chloroplast_RNA/trua4_dependent_site_overview/plot.R。
# 数据源换成当前流程：
#   WT leaf    polyARNA/STAT/fisher/{rRNA,mRNA}/WT-leaf_sites.positive.tsv
#              tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv
#   突变体     polyARNA/STAT/fisher/{rRNA,mRNA}/trua4_{1,2}-leaf_sites.Fisher.tsv
#              tRNA/STAT/fisher/ncRNA/trua4_{1,2}_Leaf_sites.Fisher.tsv
#   dependent  STAT/dependent_methylKit（is_dependent == TRUE）
# 每个 RNA 类型单独一张散点图，红点标记 trua4 dependent 位点。

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

rna_output_order <- c("rRNA", "tRNA", "mRNA")
rna_titles <- c(
  rRNA = "Chloroplast rRNA",
  tRNA = "Chloroplast tRNA",
  mRNA = "Chloroplast mRNA"
)
ratio_breaks <- seq(0, 1, by = 0.2)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]
  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }
  normalizePath(getwd())
}

read_wt_sites <- function(file_path, rna_type, chr_pattern) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "ratio")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      WT_Ratio = as.numeric(ratio)
    ) %>%
    filter(grepl(chr_pattern, chr)) %>%
    transmute(
      rna_type = factor(rna_type, levels = rna_output_order),
      chr,
      pos,
      strand,
      site_id = paste(chr, pos, strand, sep = ":"),
      WT_Ratio
    ) %>%
    arrange(rna_type, chr, pos, strand)
}

read_mutant_sites <- function(file_path, chr_pattern, ratio_col = "Treat_avg_ratio") {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", ratio_col)
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      ratio = as.numeric(.data[[ratio_col]])
    ) %>%
    filter(grepl(chr_pattern, chr)) %>%
    transmute(
      site_id = paste(chr, pos, strand, sep = ":"),
      ratio
    )
}

# 当前 methylKit 口径的 dependent-sites 表包含全部候选位点，须按 is_dependent 过滤。
read_dependent_sites <- function(file_path, chr_pattern, support_label) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "is_dependent")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      is_dependent = as.logical(is_dependent)
    ) %>%
    filter(grepl(chr_pattern, chr), !is.na(is_dependent), is_dependent) %>%
    transmute(
      site_id = paste(chr, pos, strand, sep = ":"),
      !!support_label := TRUE
    ) %>%
    distinct()
}

build_rna_table <- function(rna_type, wt_path, mut1_path, mut2_path, dep1_path, dep2_path, chr_pattern) {
  wt_tbl <- read_wt_sites(wt_path, rna_type, chr_pattern)
  mut1_tbl <- read_mutant_sites(mut1_path, chr_pattern) %>% rename(trua4_1 = ratio)
  mut2_tbl <- read_mutant_sites(mut2_path, chr_pattern) %>% rename(trua4_2 = ratio)
  dep1_tbl <- read_dependent_sites(dep1_path, chr_pattern, "dep_line1")
  dep2_tbl <- read_dependent_sites(dep2_path, chr_pattern, "dep_line2")

  wt_tbl %>%
    left_join(mut1_tbl, by = "site_id") %>%
    left_join(mut2_tbl, by = "site_id") %>%
    left_join(dep1_tbl, by = "site_id") %>%
    left_join(dep2_tbl, by = "site_id") %>%
    mutate(
      trua4_1 = as.numeric(coalesce(trua4_1, 0)),
      trua4_2 = as.numeric(coalesce(trua4_2, 0)),
      line_count = (!is.na(trua4_1) & trua4_1 > 0) + (!is.na(trua4_2) & trua4_2 > 0),
      Mut_Ratio = dplyr::case_when(
        line_count == 2 ~ (trua4_1 + trua4_2) / 2,
        line_count == 1 & trua4_1 > 0 ~ trua4_1,
        line_count == 1 & trua4_2 > 0 ~ trua4_2,
        TRUE ~ 0
      ),
      dep_line1 = coalesce(dep_line1, FALSE),
      dep_line2 = coalesce(dep_line2, FALSE),
      dependency_support = case_when(
        dep_line1 & dep_line2 ~ "shared",
        dep_line1 ~ "trua4_1_only",
        dep_line2 ~ "trua4_2_only",
        TRUE ~ "none"
      ),
      site_status = if_else(dependency_support == "none", "WT sites", "Dependent sites"),
      Delta_Ratio = WT_Ratio - Mut_Ratio
    ) %>%
    select(rna_type, chr, pos, strand, site_id, WT_Ratio, Mut_Ratio, Delta_Ratio, site_status, dependency_support)
}

build_rna_plot <- function(plot_df, rna_type, site_colors) {
  rna_df <- plot_df %>%
    filter(rna_type == !!rna_type) %>%
    mutate(site_status = factor(site_status, levels = c("WT sites", "Dependent sites")))

  if (nrow(rna_df) == 0) {
    stop("No rows found for RNA type: ", rna_type)
  }

  ggplot(rna_df, aes(x = Mut_Ratio, y = WT_Ratio)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = SIZE_TOKENS$base_linewidth,
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_point(
      data = rna_df,
      aes(fill = "WT sites"),
      size = 2.8,
      shape = 21,
      color = "black",
      stroke = 0.45
    ) +
    geom_point(
      data = filter(rna_df, site_status == "Dependent sites"),
      aes(fill = site_status),
      size = 3.0,
      shape = 21,
      color = "black",
      stroke = 0.55
    ) +
    scale_fill_manual(
      values = site_colors,
      breaks = c("WT sites", "Dependent sites"),
      limits = c("WT sites", "Dependent sites"),
      drop = FALSE,
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = ratio_breaks,
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = ratio_breaks,
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = 0.02)
    ) +
    coord_equal() +
    labs(
      title = rna_titles[[rna_type]],
      x = "trua4 Mutant Modification Ratio",
      y = "WT Modification Ratio"
    ) +
    theme(
      legend.position = c(0.98, 0.05),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      plot.title = element_text(face = "bold", hjust = 0.5, size = SIZE_TOKENS$plot_title_size),
      axis.text.x = element_text(angle = 0)
    )
}

run_analysis <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)

  site_colors <- setNames(PAL_STATUS, c("WT sites", "Dependent sites"))

  specs <- list(
    list(
      rna_type = "rRNA",
      wt_path = file.path(repo_root, "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv"),
      mut1_path = file.path(repo_root, "polyARNA/STAT/fisher/rRNA/trua4_1-leaf_sites.Fisher.tsv"),
      mut2_path = file.path(repo_root, "polyARNA/STAT/fisher/rRNA/trua4_2-leaf_sites.Fisher.tsv"),
      dep1_path = file.path(repo_root, "polyARNA/STAT/dependent_methylKit/rRNA/leaf/trua4/trua4_1-leaf.dependent-sites.tsv"),
      dep2_path = file.path(repo_root, "polyARNA/STAT/dependent_methylKit/rRNA/leaf/trua4/trua4_2-leaf.dependent-sites.tsv"),
      chr_pattern = "^Ch23S$"
    ),
    list(
      rna_type = "tRNA",
      wt_path = file.path(repo_root, "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv"),
      mut1_path = file.path(repo_root, "tRNA/STAT/fisher/ncRNA/trua4_1_Leaf_sites.Fisher.tsv"),
      mut2_path = file.path(repo_root, "tRNA/STAT/fisher/ncRNA/trua4_2_Leaf_sites.Fisher.tsv"),
      dep1_path = file.path(repo_root, "tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_1-leaf.dependent-sites.tsv"),
      dep2_path = file.path(repo_root, "tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_2-leaf.dependent-sites.tsv"),
      chr_pattern = "^tRNA-Pt-"
    ),
    list(
      rna_type = "mRNA",
      wt_path = file.path(repo_root, "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv"),
      mut1_path = file.path(repo_root, "polyARNA/STAT/fisher/mRNA/trua4_1-leaf_sites.Fisher.tsv"),
      mut2_path = file.path(repo_root, "polyARNA/STAT/fisher/mRNA/trua4_2-leaf_sites.Fisher.tsv"),
      dep1_path = file.path(repo_root, "polyARNA/STAT/dependent_methylKit/mRNA/leaf/trua4/trua4_1-leaf.dependent-sites.tsv"),
      dep2_path = file.path(repo_root, "polyARNA/STAT/dependent_methylKit/mRNA/leaf/trua4/trua4_2-leaf.dependent-sites.tsv"),
      chr_pattern = "^ATCG"
    )
  )

  missing_inputs <- unique(unlist(lapply(specs, function(spec) {
    paths <- unname(unlist(spec[c("wt_path", "mut1_path", "mut2_path", "dep1_path", "dep2_path")]))
    paths[!file.exists(paths)]
  })))
  if (length(missing_inputs) > 0) {
    stop("Missing required inputs:\n", paste(missing_inputs, collapse = "\n"))
  }

  plot_tbl <- bind_rows(lapply(specs, function(spec) {
    build_rna_table(
      rna_type = spec$rna_type,
      wt_path = spec$wt_path,
      mut1_path = spec$mut1_path,
      mut2_path = spec$mut2_path,
      dep1_path = spec$dep1_path,
      dep2_path = spec$dep2_path,
      chr_pattern = spec$chr_pattern
    )
  })) %>%
    mutate(
      rna_type = factor(rna_type, levels = rna_output_order),
      site_status = factor(site_status, levels = c("WT sites", "Dependent sites"))
    ) %>%
    arrange(rna_type, chr, pos, strand)

  if (nrow(plot_tbl) == 0) {
    stop("No chloroplast WT sites were available for the trua4-dependent overview.")
  }

  summary_tbl <- plot_tbl %>%
    group_by(rna_type) %>%
    summarise(
      n_sites = n(),
      n_dependent = sum(site_status == "Dependent sites"),
      label = paste0("n = ", scales::comma(n_sites), "; Dependent = ", scales::comma(n_dependent)),
      .groups = "drop"
    )

  readr::write_tsv(plot_tbl, file.path(script_dir, "trua4_dependent_site_overview_dotplot.tsv"))
  readr::write_tsv(summary_tbl, file.path(script_dir, "trua4_dependent_site_overview_summary.tsv"))

  for (rna_name in rna_output_order) {
    stem <- file.path(script_dir, paste0("trua4_dependent_site_overview_dotplot_", tolower(rna_name)))
    save_plot(build_rna_plot(plot_tbl, rna_name, site_colors), stem, width = 8, height = 6.6)
  }
}

if (sys.nframe() == 0) {
  run_analysis()
}
