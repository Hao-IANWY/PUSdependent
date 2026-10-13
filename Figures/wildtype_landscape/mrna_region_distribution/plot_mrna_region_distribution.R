#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggprism)
  library(sysfonts)
  library(showtext)
  library(tibble)
})

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

normalize_region <- function(annotation) {
  dplyr::case_when(
    annotation == "five_prime_UTR" ~ "5'UTR",
    annotation == "three_prime_UTR" ~ "3'UTR",
    annotation == "CDS" ~ "CDS",
    TRUE ~ "Other"
  )
}

read_annotated_sites <- function(file_path, tissue_name) {
  readr::read_tsv(file_path, show_col_types = FALSE) |>
    transmute(
      Tissue = tissue_name,
      Region = normalize_region(annoation)
    )
}

run_fig1d <- function(out_dir = NULL, output_stem = "mrna_region_distribution", save_outputs = TRUE, fig1d_dir = NULL) {
  script_dir <- get_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "../../.."), mustWork = TRUE)
  source(file.path(repo_root, "Figures/shared/plot_style.R"))
  source(file.path(repo_root, "Figures/shared/plot_tokens.R"))
  init_theme(base_size = 16, border = TRUE, sizes = SIZE_TOKENS)
  if (is.null(fig1d_dir)) {
    fig1d_dir <- script_dir
  }
  site_dir <- file.path(fig1d_dir, "site_analyse")
  out_dir <- if (is.null(out_dir)) script_dir else normalizePath(out_dir, mustWork = FALSE)

  # positive: Ψ 阳性位点(分子); expressed_pool: 表达(覆盖达标)的 U 位点全集(分母)。
  # 注意 STAT/fisher 的 Fisher.tsv 已被 Pseudo_callsites.R 按 pooled p.adj < 0.01 过滤,
  # 并非全部候选位点; 表达池由 build_expressed_pool.sh 按同一入池规则从 filter_T.tsv 重建
  # (每重复 depth >= 10 的 T 位点 6 重复并集, 缺失记 0, Mock/Treat 组均 depth >= 10)。
  # 分母即该组织"表达了"的区域序列量, 校正了样本间测序深度/降解差异(read bias),
  # 回应审稿人对 Fig 1f 的意见。
  fisher_dir <- file.path(repo_root, "polyARNA/STAT/fisher/mRNA")
  input_files <- tibble::tribble(
    ~Tissue,  ~prefix,
    "Root",   "WT-Root",
    "Leaf",   "WT-leaf",
    "Flower", "WT-Flower"
  ) |>
    mutate(
      positive_file = file.path(site_dir, paste0(prefix, "_sites.positive.annoated.tsv")),
      expressed_file = file.path(site_dir, paste0(prefix, ".expressed_pool.annoated.tsv")),
      fisher_file = file.path(fisher_dir, paste0(prefix, "_sites.Fisher.tsv"))
    )

  missing_files <- c(input_files$positive_file, input_files$expressed_file, input_files$fisher_file)
  missing_files <- missing_files[!file.exists(missing_files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing annotated input files:\n",
      paste(missing_files, collapse = "\n"),
      "\nRun annotate_mrna_positive_sites.sh and build_expressed_pool.sh first."
    )
  }

  display_regions <- c("5'UTR", "CDS", "3'UTR")
  tissue_levels <- c("Root", "Leaf", "Flower")

  positive_counts <- purrr::map2_dfr(input_files$positive_file, input_files$Tissue, read_annotated_sites) |>
    count(Tissue, Region, name = "Count")
  expressed_counts <- purrr::map2_dfr(input_files$expressed_file, input_files$Tissue, read_annotated_sites) |>
    count(Tissue, Region, name = "Expressed")

  summary_df <- tidyr::expand_grid(Tissue = tissue_levels, Region = display_regions) |>
    left_join(positive_counts, by = c("Tissue", "Region")) |>
    left_join(expressed_counts, by = c("Tissue", "Region")) |>
    mutate(
      Count = tidyr::replace_na(Count, 0L),
      Expressed = tidyr::replace_na(Expressed, 0L)
    )

  if (any(summary_df$Expressed == 0)) {
    stop("At least one tissue x region combination has zero expressed positions; cannot compute rates.")
  }

  summary_df <- summary_df |>
    mutate(
      Rate_per_1000 = Count / Expressed * 1000,
      Tissue = factor(Tissue, levels = tissue_levels),
      Region = factor(Region, levels = display_regions)
    )

  # 2x3 列联检验(阳性/未阳性 x 组织): 优先 Fisher 精确检验;
  # 计数过大超出精确检验工作区(FEXACT error 40)时回退卡方, 此时期望计数足够大, 近似可靠
  test_rate_across_tissues <- function(count, total) {
    tab <- matrix(c(count, total - count), nrow = 3)
    tryCatch(
      list(method = "Fisher's exact test", p_value = stats::fisher.test(tab)$p.value),
      error = function(e) list(method = "Chi-squared test", p_value = stats::chisq.test(tab)$p.value)
    )
  }

  # 按区域检验三个组织的修饰率是否一致
  stats_df <- purrr::map_dfr(display_regions, function(region) {
    region_df <- dplyr::filter(summary_df, Region == region)
    test <- test_rate_across_tissues(region_df$Count, region_df$Expressed)
    tibble::tibble(Region = region, method = test$method, p_value = test$p_value)
  }) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"))

  # ---- 每 1000 nt 注释区域长度的修饰密度(审稿人字面口径, 与 per-U 版并列) ----
  # 分母为 Araport11 代表转录本集合的区域拼接总长(常数, 不随组织变),
  # 由 compute_region_lengths.sh 生成
  region_lengths <- readr::read_tsv(
    file.path(script_dir, "mrna_region_distribution_region_lengths.tsv"),
    show_col_types = FALSE
  )
  density_df <- summary_df |>
    transmute(Tissue, Region = as.character(Region), Count) |>
    left_join(region_lengths, by = "Region") |>
    mutate(
      Density_per_1000nt = Count / Length_nt * 1000,
      Region = factor(Region, levels = display_regions)
    )

  density_stats_df <- purrr::map_dfr(display_regions, function(region) {
    region_df <- dplyr::filter(density_df, Region == region)
    test <- test_rate_across_tissues(region_df$Count, region_df$Length_nt)
    tibble::tibble(Region = region, method = test$method, p_value = test$p_value)
  }) |>
    mutate(p_adj = p.adjust(p_value, method = "BH"))

  # 按重复的分区域修饰率(供 response letter: 正面回应审稿人对样本间差异的关切)。
  # 分子 = RepX_pass 标志位点数(注意: 这些标志只在 pooled 显著的候选位点上计算);
  # 分母与主图一致(表达池)。重复数随流程配置变化(3 -> 4), pass 列按表头动态识别,
  # 不硬编码 Rep1-3。
  per_replicate_long <- purrr::map_dfr(seq_len(nrow(input_files)), function(i) {
    sites <- readr::read_tsv(input_files$expressed_file[i], show_col_types = FALSE) |>
      transmute(chr, pos, Region = normalize_region(annoation))
    pass_cols <- grep("^Rep\\d+_pass$", names(readr::read_tsv(
      input_files$fisher_file[i], show_col_types = FALSE, n_max = 0)), value = TRUE)
    readr::read_tsv(input_files$fisher_file[i], show_col_types = FALSE,
                    col_select = c("chr", "pos", dplyr::all_of(pass_cols))) |>
      inner_join(sites, by = c("chr", "pos")) |>
      filter(Region %in% display_regions) |>
      tidyr::pivot_longer(dplyr::all_of(pass_cols), names_to = "Replicate", values_to = "Pass") |>
      mutate(Replicate = sub("_pass$", "", Replicate)) |>
      count(Region, Replicate, wt = Pass == 1, name = "Count") |>
      mutate(Tissue = input_files$Tissue[i])
  }) |>
    left_join(
      summary_df |>
        transmute(Tissue = as.character(Tissue), Region = as.character(Region), Expressed),
      by = c("Tissue", "Region")
    ) |>
    mutate(Rate_per_1000 = Count / Expressed * 1000) |>
    relocate(Tissue, Region, Replicate, Count, Expressed, Rate_per_1000) |>
    arrange(match(Tissue, tissue_levels), match(Region, display_regions), Replicate)

  per_replicate_wide <- per_replicate_long |>
    select(Tissue, Region, Replicate, Rate_per_1000) |>
    tidyr::pivot_wider(names_from = Replicate, values_from = Rate_per_1000, names_glue = "{Replicate}_rate") |>
    rowwise() |>
    mutate(
      Mean_rate = mean(c_across(ends_with("_rate"))),
      SD_rate = sd(c_across(ends_with("_rate")))
    ) |>
    ungroup() |>
    arrange(match(Tissue, tissue_levels), match(Region, display_regions))

  fill_values <- c(
    "5'UTR" = "#4DBBD5CC",
    "CDS"   = "#E64B35CC",
    "3'UTR" = "#3C5488CC"
  )

  dodge <- position_dodge(width = 0.78)

  # per-U 修饰率图与 per-nt 密度图共用同一样式
  build_region_barplot <- function(plot_df, y_column, y_label) {
    ggplot(plot_df, aes(x = Tissue, y = .data[[y_column]], fill = Region)) +
      geom_col(
        width = 0.72,
        color = "black",
        linewidth = 0.25,
        position = dodge
      ) +
      scale_fill_manual(
        values = fill_values,
        breaks = display_regions,
        name = NULL,
        drop = TRUE
      ) +
      scale_y_continuous(
        expand = expansion(mult = c(0, 0.08)),
        guide = ggprism::guide_prism_offset()
      ) +
      # 内部 Tissue 保持单数键(下游 join 依赖), 仅 x 轴显示层改复数
      scale_x_discrete(name = NULL, drop = TRUE,
                       labels = c(Root = "Roots", Leaf = "Leaves", Flower = "Flowers")) +
      labs(
        x = NULL,
        y = y_label
      ) +
      coord_cartesian(clip = "off") +
      theme(
        text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = SIZE_TOKENS$base_linewidth),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = margin(18, 18, 14, 18)
      )
  }

  bar_plot <- build_region_barplot(summary_df, "Rate_per_1000", "Ψ sites per 1,000 expressed uridines")
  density_plot <- build_region_barplot(density_df, "Density_per_1000nt", "Ψ sites per 1,000 nt")

  outputs <- list(
    summary_tsv = file.path(out_dir, paste0(output_stem, "_rates.tsv")),
    stats_tsv = file.path(out_dir, paste0(output_stem, "_stats.tsv")),
    density_tsv = file.path(out_dir, paste0(output_stem, "_density.tsv")),
    density_stats_tsv = file.path(out_dir, paste0(output_stem, "_density_stats.tsv")),
    per_replicate_tsv = file.path(out_dir, paste0(output_stem, "_per_replicate_rates.tsv")),
    per_replicate_summary_tsv = file.path(out_dir, paste0(output_stem, "_per_replicate_summary.tsv")),
    plot_pdf = file.path(out_dir, paste0(output_stem, "_barplot.pdf")),
    plot_svg = file.path(out_dir, paste0(output_stem, "_barplot.svg")),
    plot_png = file.path(out_dir, paste0(output_stem, "_barplot.png")),
    density_pdf = file.path(out_dir, paste0(output_stem, "_density_barplot.pdf")),
    density_svg = file.path(out_dir, paste0(output_stem, "_density_barplot.svg")),
    density_png = file.path(out_dir, paste0(output_stem, "_density_barplot.png"))
  )

  if (save_outputs) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(summary_df, outputs$summary_tsv)
    readr::write_tsv(stats_df, outputs$stats_tsv)
    readr::write_tsv(density_df, outputs$density_tsv)
    readr::write_tsv(density_stats_df, outputs$density_stats_tsv)
    readr::write_tsv(per_replicate_long, outputs$per_replicate_tsv)
    readr::write_tsv(per_replicate_wide, outputs$per_replicate_summary_tsv)
    save_plot(bar_plot, file.path(out_dir, paste0(output_stem, "_barplot")), width = 5.5, height = 6.5)
    save_plot(density_plot, file.path(out_dir, paste0(output_stem, "_density_barplot")), width = 5.5, height = 6.5)

    message("Saved summary table to: ", outputs$summary_tsv)
    message("Saved region stats to: ", outputs$stats_tsv)
    message("Saved density table to: ", outputs$density_tsv)
    message("Saved per-replicate rates to: ", outputs$per_replicate_tsv)
    message("Saved plot files to: ", out_dir)
  }

  list(
    plot = bar_plot,
    density_plot = density_plot,
    summary_table = summary_df,
    stats_table = stats_df,
    density_table = density_df,
    density_stats_table = density_stats_df,
    per_replicate_rates = per_replicate_long,
    per_replicate_summary = per_replicate_wide,
    output_paths = outputs,
    out_dir = out_dir,
    output_stem = output_stem
  )
}

if (sys.nframe() == 0) {
  invisible(run_fig1d())
}
