#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(DOSE)
  library(ggplot2)
  library(tidyr)
  library(forcats)
  library(ggprism)
  library(stringr)
  library(org.At.tair.db)
  library(sysfonts)
  library(showtext)
  library(patchwork)
})

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  script_path <- cmd_args[grepl(file_arg, cmd_args)]
  if (length(script_path) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", script_path[1]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
setwd(script_dir)

org_db <- org.At.tair.db
export_dpi <- 600
pdf_font_family <- "Arial"
arial_font_path <- "/home/HL/miniforge3/envs/analysis/fonts/arial.ttf"
if (!file.exists(arial_font_path)) {
  stop(
    "Default Arial font file not found at ", arial_font_path, ". ",
    "Please locate an Arial font file with: fc-list | grep -i 'arial'"
  )
}

sysfonts::font_add(pdf_font_family, regular = arial_font_path)
showtext::showtext_opts(dpi = export_dpi)

theme_set(
  ggprism::theme_prism(base_size = 24, border = TRUE) +
    theme(
      text = element_text(family = pdf_font_family, face = "bold", colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.6),
      legend.key.spacing.y = grid::unit(12, "points"),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      strip.background = element_blank(),
      strip.placement = "outside"
    )
)

gene_files <- c(
  Root_unique = file.path(script_dir, "root_unique_gene_list.tsv"),
  Leaf_unique = file.path(script_dir, "leaf_unique_gene_list.tsv"),
  Flower_unique = file.path(script_dir, "flower_unique_gene_list.tsv"),
  Common = file.path(script_dir, "common_gene_list.tsv")
)

read_gene_list <- function(gene_file) {
  if (!file.exists(gene_file)) {
    stop("Missing gene list file: ", gene_file)
  }
  genes <- read.table(
    gene_file,
    header = FALSE,
    sep = "\t",
    stringsAsFactors = FALSE,
    colClasses = "character"
  )
  if (nrow(genes) == 0) return(character(0))
  sort(unique(genes[[1]]))
}

write_empty_go_table <- function(prefix) {
  write.table(
    data.frame(),
    file = file.path(script_dir, paste0(prefix, "_GO.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

collect_go_terms <- function(gene_ids, prefix) {
  result <- list(
    mapped_genes = 0L,
    significant_terms = 0L,
    plotted_terms = 0L,
    plot_df = NULL
  )

  if (length(gene_ids) == 0) {
    write_empty_go_table(prefix)
    return(result)
  }

  df_entrez <- clusterProfiler::bitr(
    gene_ids,
    fromType = "TAIR",
    toType = "ENTREZID",
    OrgDb = org_db
  )
  result$mapped_genes <- nrow(distinct(df_entrez, ENTREZID, .keep_all = TRUE))

  if (nrow(df_entrez) == 0) {
    write_empty_go_table(prefix)
    return(result)
  }

  df_entrez <- distinct(df_entrez, ENTREZID, .keep_all = TRUE)

  go_enrich <- clusterProfiler::enrichGO(
    gene = df_entrez$ENTREZID,
    ont = "all",
    keyType = "ENTREZID",
    OrgDb = org_db,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05
  )

  if (is.null(go_enrich) || nrow(as.data.frame(go_enrich)) == 0) {
    write_empty_go_table(prefix)
    return(result)
  }

  go_enrich <- DOSE::setReadable(
    go_enrich,
    OrgDb = org_db,
    keyType = "ENTREZID"
  )

  go_simplified <- clusterProfiler::simplify(
    go_enrich,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min
  )

  if (is.null(go_simplified) || nrow(as.data.frame(go_simplified)) == 0) {
    write_empty_go_table(prefix)
    return(result)
  }

  go_result <- go_simplified@result
  result$significant_terms <- nrow(go_result)

  if (nrow(go_result) == 0) {
    write.table(
      go_result,
      file = file.path(script_dir, paste0(prefix, "_GO.tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    return(result)
  }

  go_result <- go_result |>
    tidyr::separate(
      col = GeneRatio,
      into = c("GR1", "GR2"),
      sep = "/",
      remove = FALSE
    ) |>
    dplyr::mutate(GeneRatio = as.numeric(GR1) / as.numeric(GR2))

  write.table(
    go_result,
    file = file.path(script_dir, paste0(prefix, "_GO.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  plot_df <- dplyr::bind_rows(
    utils::head(dplyr::filter(go_result, ONTOLOGY == "BP"), 5L),
    utils::head(dplyr::filter(go_result, ONTOLOGY == "CC"), 5L),
    utils::head(dplyr::filter(go_result, ONTOLOGY == "MF"), 5L)
  )

  if (nrow(plot_df) > 0) {
    plot_df <- plot_df |>
      dplyr::mutate(
        ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF")),
        Description = stringr::str_wrap(Description, width = 28)
      )
    result$plot_df <- plot_df
    result$plotted_terms <- nrow(plot_df)
  }

  result
}

build_single_plot <- function(plot_df, panel_title, fill_limits, size_limits, size_breaks) {
  ggplot(plot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio))) +
    geom_point(
      aes(size = Count, fill = p.adjust),
      shape = 21,
      color = "black",
      stroke = 0.6,
      alpha = 0.8
    ) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.10))) +
    scale_fill_gradient(
      low = "#E64B35CC",
      high = "#3C5488CC",
      limits = fill_limits,
      guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(
      name = "Count",
      range = c(4, 10),
      limits = size_limits,
      breaks = size_breaks
    ) +
    labs(
      title = panel_title,
      x = "Gene ratio",
      y = NULL,
      fill = "Adjusted p",
      size = "Count"
    ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b = 10, unit = "points")),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_text(size = 24),
      strip.text.y = element_text(size = 20),
      panel.spacing.y = grid::unit(0.8, "lines"),
      legend.position = "right",
      legend.title = element_text(size = 18),
      legend.text = element_text(size = 16),
      legend.key.height = grid::unit(1.2, "lines"),
      legend.key.width = grid::unit(1.2, "lines")
    )
}

build_empty_plot <- function(panel_title) {
  ggplot() +
    annotate(
      "text",
      x = 0.5,
      y = 0.5,
      label = "No enriched GO terms\n(BH-adjusted P < 0.05, q < 0.05)",
      family = pdf_font_family,
      fontface = "bold",
      size = 6.5,
      color = "#4D4D4D"
    ) +
    labs(title = panel_title) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      text = element_text(family = pdf_font_family, face = "bold", colour = "black"),
      plot.title = element_text(size = 20, hjust = 0.5, margin = margin(b = 10, unit = "points")),
      plot.margin = margin(18, 18, 18, 18)
    )
}

save_single_plot <- function(plot_obj, prefix) {
  ggplot2::ggsave(
    filename = file.path(script_dir, paste0(prefix, "_go_enrichment.pdf")),
    plot = plot_obj,
    width = 8,
    height = 10,
    dpi = export_dpi,
    device = function(filename, width, height, ...) {
      grDevices::cairo_pdf(
        filename = filename,
        width = width,
        height = height,
        family = pdf_font_family,
        ...
      )
    }
  )
  ggplot2::ggsave(
    filename = file.path(script_dir, paste0(prefix, "_go_enrichment.png")),
    plot = plot_obj,
    width = 8,
    height = 10,
    dpi = 220,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(script_dir, paste0(prefix, "_go_enrichment.svg")),
    plot = plot_obj,
    width = 8,
    height = 10,
    device = grDevices::svg,
    dpi = export_dpi
  )
}

gene_lists <- lapply(gene_files, read_gene_list)
group_name_map <- vapply(
  names(gene_lists),
  function(prefix) {
    sprintf("%s ( n = %d )", gsub("_", "-", prefix), length(gene_lists[[prefix]]))
  },
  character(1)
)

analysis_results <- lapply(names(gene_lists), function(prefix) {
  collect_go_terms(gene_lists[[prefix]], prefix)
})
names(analysis_results) <- names(gene_lists)

all_plot_data <- bind_rows(lapply(analysis_results, `[[`, "plot_df"))
if (nrow(all_plot_data) > 0) {
  fill_limits <- range(all_plot_data$p.adjust, na.rm = TRUE)
  size_limits <- range(all_plot_data$Count, na.rm = TRUE)
  if (diff(fill_limits) == 0) fill_limits <- fill_limits + c(-0.001, 0.001)
  if (diff(size_limits) == 0) size_limits <- size_limits + c(-0.5, 0.5)
  size_breaks <- unique(pretty(size_limits, n = 4))
  size_breaks <- size_breaks[size_breaks >= size_limits[1] & size_breaks <= size_limits[2]]
  if (length(size_breaks) == 0) size_breaks <- round(size_limits)
} else {
  fill_limits <- c(0, 0.05)
  size_limits <- c(0, 1)
  size_breaks <- c(0, 1)
}

plot_objects <- lapply(names(analysis_results), function(prefix) {
  result <- analysis_results[[prefix]]
  if (is.null(result$plot_df) || nrow(result$plot_df) == 0) {
    plot_obj <- build_empty_plot(group_name_map[[prefix]])
  } else {
    plot_obj <- build_single_plot(
      result$plot_df,
      group_name_map[[prefix]],
      fill_limits,
      size_limits,
      size_breaks
    )
  }
  save_single_plot(plot_obj, prefix)
  plot_obj
})
names(plot_objects) <- names(analysis_results)

combined_plot <- (
  patchwork::wrap_plots(plotlist = plot_objects, nrow = 1) |
    patchwork::guide_area()
) +
  patchwork::plot_layout(
    widths = c(rep(1, length(plot_objects)), 0.32),
    guides = "collect"
  ) &
  theme(legend.position = "right")

combined_outputs <- list(
  pdf = file.path(script_dir, "mrna_go_enrichment_combined.pdf"),
  png = file.path(script_dir, "mrna_go_enrichment_combined.png"),
  svg = file.path(script_dir, "mrna_go_enrichment_combined.svg")
)

ggplot2::ggsave(
  filename = combined_outputs$pdf,
  plot = combined_plot,
  width = 32,
  height = 10,
  dpi = export_dpi,
  limitsize = FALSE,
  device = function(filename, width, height, ...) {
    grDevices::cairo_pdf(
      filename = filename,
      width = width,
      height = height,
      family = pdf_font_family,
      ...
    )
  }
)
ggplot2::ggsave(
  filename = combined_outputs$png,
  plot = combined_plot,
  width = 32,
  height = 10,
  dpi = 220,
  limitsize = FALSE,
  bg = "white"
)
ggplot2::ggsave(
  filename = combined_outputs$svg,
  plot = combined_plot,
  width = 32,
  height = 10,
  dpi = export_dpi,
  limitsize = FALSE,
  device = grDevices::svg
)

go_summary <- dplyr::tibble(
  Group = names(analysis_results),
  Input_genes = vapply(gene_lists, length, integer(1)),
  Mapped_Entrez_genes = vapply(analysis_results, function(x) x$mapped_genes, integer(1)),
  Significant_GO_terms = vapply(analysis_results, function(x) x$significant_terms, integer(1)),
  Plotted_GO_terms = vapply(analysis_results, function(x) x$plotted_terms, integer(1))
)
write.table(
  go_summary,
  file = file.path(script_dir, "mrna_go_enrichment_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("[GO enrichment] Analysis complete.\n")
print(as.data.frame(go_summary), row.names = FALSE)
