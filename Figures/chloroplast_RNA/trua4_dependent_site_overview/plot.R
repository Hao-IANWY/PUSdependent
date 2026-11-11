#!/usr/bin/env Rscript
# trua4_dependent_site_overview — WT vs trua4 scatter plots of chloroplast
# rRNA / tRNA / mRNA sites in leaf
#
# Inputs (beta-binomial calling criteria, STAT/dependent_betabinomial):
#   WT sites and modification fractions   <rna>/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv
#   Per-allele fractions and calls        <rna>/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv
#
# The dependency calling criteria match Figures/dependent_betabinomial: a
# site is highlighted as dependent only when is_dependent is TRUE in both
# mutant alleles; sites supported by a single allele are recorded only in
# dependency_support and are not highlighted.
#
# Chloroplast filters: rRNA chr == Ch23S; tRNA chr ^tRNA-Pt-; mRNA chr ^ATCG.
# Each RNA class gets one scatter plot (WT fraction vs mutant fraction) and a
# summary TSV.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

ratio_breaks <- seq(0, 1, by = 0.2)
rna_output_order <- c("rRNA", "tRNA", "mRNA")
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

read_wt_sites <- function(file_path, rna_type, chr_pattern) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "Treat_avg_ratio")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand),
      WT_Ratio = as.numeric(Treat_avg_ratio)
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

read_line_sites <- function(file_path, line_name, chr_pattern) {
  dat <- readr::read_tsv(file_path, show_col_types = FALSE)
  required_cols <- c("chr", "pos", "strand", "Mut_Ratio", "is_dependent")
  missing_cols <- setdiff(required_cols, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", file_path, ": ", paste(missing_cols, collapse = ", "))
  }

  dat %>%
    mutate(
      chr = as.character(chr),
      pos = as.integer(pos),
      strand = as.character(strand)
    ) %>%
    filter(grepl(chr_pattern, chr)) %>%
    transmute(
      site_id = paste(chr, pos, strand, sep = ":"),
      !!line_name := as.numeric(Mut_Ratio),
      !!paste0(line_name, "_dep") := as.logical(is_dependent)
    )
}

build_rna_table <- function(spec) {
  wt_tbl <- read_wt_sites(spec$wt_path, spec$rna_type, spec$chr_pattern)
  line1_tbl <- read_line_sites(spec$mut1_path, "trua4_1", spec$chr_pattern)
  line2_tbl <- read_line_sites(spec$mut2_path, "trua4_2", spec$chr_pattern)

  wt_tbl %>%
    left_join(line1_tbl, by = "site_id") %>%
    left_join(line2_tbl, by = "site_id") %>%
    mutate(
      trua4_1 = as.numeric(trua4_1),
      trua4_2 = as.numeric(trua4_2),
      trua4_1_dep = coalesce(trua4_1_dep, FALSE),
      trua4_2_dep = coalesce(trua4_2_dep, FALSE),
      Mut_Ratio = rowMeans(cbind(trua4_1, trua4_2), na.rm = TRUE),
      Mut_Ratio = if_else(is.nan(Mut_Ratio), 0, Mut_Ratio),
      dependency_support = case_when(
        trua4_1_dep & trua4_2_dep ~ "shared",
        trua4_1_dep ~ "trua4_1_only",
        trua4_2_dep ~ "trua4_2_only",
        TRUE ~ "none"
      ),
      site_status = if_else(
        trua4_1_dep & trua4_2_dep,
        "Dependent sites",
        "WT sites"
      ),
      Delta_Ratio = WT_Ratio - Mut_Ratio
    ) %>%
    select(
      rna_type, chr, pos, strand, site_id,
      WT_Ratio, Mut_Ratio, Delta_Ratio, site_status, dependency_support,
      trua4_1, trua4_2, trua4_1_dep, trua4_2_dep
    )
}

build_rna_plot <- function(plot_df, rna_name) {
  rna_df <- plot_df %>%
    filter(rna_type == !!rna_name) %>%
    mutate(site_status = factor(site_status, levels = c("WT sites", "Dependent sites")))

  if (nrow(rna_df) == 0) {
    stop("No rows found for RNA type: ", rna_name)
  }

  status_colors <- setNames(PAL_STATUS, c("WT sites", "Dependent sites"))

  ggplot(rna_df, aes(x = Mut_Ratio, y = WT_Ratio)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linewidth = SIZE_TOKENS$base_linewidth,
      linetype = "dashed",
      color = "grey45"
    ) +
    geom_point(
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
      values = status_colors,
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
      title = rna_titles[[rna_name]],
      x = "*trua4* Mutant Modification Ratio",
      y = "WT Modification Ratio"
    ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.direction = "vertical",
      # theme_prism hard-codes legend.text to plain; restore Arial Bold
      legend.text = element_text(
        family = "Arial",
        face = "bold",
        colour = "black",
        size = SIZE_TOKENS$legend_text_size
      ),
      plot.title = element_text(hjust = 0.5, size = SIZE_TOKENS$plot_title_size),
      axis.title.x = ggtext::element_markdown(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      plot.margin = margin(10, 20, 10, 12)
    )
}

run_analysis <- function() {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = SIZE_TOKENS$base_size, border = TRUE, sizes = SIZE_TOKENS)

  polyarna_leaf <- file.path(repo_root, "polyARNA/STAT/dependent_betabinomial")
  trna_leaf <- file.path(repo_root, "tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4")

  specs <- list(
    list(
      rna_type = "rRNA",
      wt_path = file.path(polyarna_leaf, "rRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv"),
      mut1_path = file.path(polyarna_leaf, "rRNA/leaf/trua4/trua4_1-leaf.dependent-sites.tsv"),
      mut2_path = file.path(polyarna_leaf, "rRNA/leaf/trua4/trua4_2-leaf.dependent-sites.tsv"),
      chr_pattern = "^Ch23S$"
    ),
    list(
      rna_type = "tRNA",
      wt_path = file.path(trna_leaf, "WT-leaf_sites.Treat_avg_ratio.tsv"),
      mut1_path = file.path(trna_leaf, "trua4_1-leaf.dependent-sites.tsv"),
      mut2_path = file.path(trna_leaf, "trua4_2-leaf.dependent-sites.tsv"),
      chr_pattern = "^tRNA-Pt-"
    ),
    list(
      rna_type = "mRNA",
      wt_path = file.path(polyarna_leaf, "mRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv"),
      mut1_path = file.path(polyarna_leaf, "mRNA/leaf/trua4/trua4_1-leaf.dependent-sites.tsv"),
      mut2_path = file.path(polyarna_leaf, "mRNA/leaf/trua4/trua4_2-leaf.dependent-sites.tsv"),
      chr_pattern = "^ATCG"
    )
  )

  missing_inputs <- unique(unlist(lapply(specs, function(spec) {
    paths <- unname(unlist(spec[c("wt_path", "mut1_path", "mut2_path")]))
    paths[!file.exists(paths)]
  })))
  if (length(missing_inputs) > 0) {
    stop("Missing required inputs:\n", paste(missing_inputs, collapse = "\n"))
  }

  plot_tbl <- bind_rows(lapply(specs, build_rna_table)) %>%
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
      label = paste0("n = ", n_sites, "; Dependent = ", n_dependent),
      .groups = "drop"
    )

  readr::write_tsv(
    plot_tbl,
    file.path(script_dir, "trua4_dependent_site_overview_dotplot.tsv")
  )
  readr::write_tsv(
    summary_tbl,
    file.path(script_dir, "trua4_dependent_site_overview_summary.tsv")
  )

  for (rna_name in rna_output_order) {
    stem <- file.path(
      script_dir,
      paste0("trua4_dependent_site_overview_dotplot_", tolower(rna_name))
    )
    save_plot(
      build_rna_plot(plot_tbl, rna_name),
      stem,
      width = 8,
      height = 6.6
    )
  }
}

if (sys.nframe() == 0) {
  run_analysis()
}
