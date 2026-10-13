#!/usr/bin/env Rscript
# plot.R — Region proportions of TRUA-dependent mRNA sites by tissue and genotype.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))
init_theme()  # 原依赖 theme.R source 时的顶层 theme_set 副作用(base_size=16, border=TRUE)
base_linewidth <- SIZE_TOKENS$base_linewidth

normalize_region <- function(annotation) {
  dplyr::case_when(
    annotation == "five_prime_UTR" ~ "5'UTR",
    annotation == "three_prime_UTR" ~ "3'UTR",
    annotation == "CDS" ~ "CDS",
    TRUE ~ "Other"
  )
}

parse_file_metadata <- function(file_path) {
  base_name <- basename(file_path)
  matches <- regexec("^([^-]+)-(.+)\\.dependent\\.annoated\\.tsv$", base_name)
  parts <- regmatches(base_name, matches)[[1]]
  if (length(parts) != 3) {
    stop("Cannot parse genotype and tissue from file name: ", base_name)
  }
  list(genotype = parts[2], tissue = parts[3])
}

read_annotated_sites <- function(file_path) {
  metadata <- parse_file_metadata(file_path)
  tbl <- read.table(
    file_path,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!"annoation" %in% colnames(tbl)) {
    stop("Annotation column 'annoation' not found in: ", file_path)
  }
  # methylKit can produce header-only files (zero dependent sites) for some
  # genotype × tissue panels. Return an empty frame so they survive the
  # downstream rbind and are dropped later at the proportion step.
  if (nrow(tbl) == 0) {
    return(data.frame(
      Genotype = character(0),
      Tissue = character(0),
      Region = character(0),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    Genotype = metadata$genotype,
    Tissue = metadata$tissue,
    Region = normalize_region(tbl$annoation),
    stringsAsFactors = FALSE
  )
}

annotated_files <- list.files(
  file.path(script_dir, "site_analyse"),
  pattern = "^trua.*-.+\\.dependent\\.annoated\\.tsv$",
  full.names = TRUE
)
if (length(annotated_files) != 21) {
  stop("Expected 21 annotated dependent files, found ", length(annotated_files), ". Run data_prep.sh first.")
}

plot_df <- do.call(rbind, lapply(annotated_files, read_annotated_sites)) |>
  mutate(
    Tissue = factor(Tissue, levels = c("Root", "leaf", "Flower"),
                    labels = c("Root", "Leaf", "Flower")),
    Genotype = factor(Genotype, levels = c("trua1a", "trua1b", "trua2",
                                           "trua3", "trua4", "trua5", "trua6")),
    Region = factor(Region, levels = c("5'UTR", "CDS", "3'UTR", "Other"))
  )

summary_df <- plot_df |>
  count(Tissue, Genotype, Region, name = "Count") |>
  complete(Tissue, Genotype, Region, fill = list(Count = 0)) |>
  filter(as.character(Region) %in% c("5'UTR", "CDS", "3'UTR")) |>
  group_by(Tissue, Genotype) |>
  mutate(
    Total = sum(Count),
    Proportion = Count / Total
  ) |>
  ungroup() |>
  mutate(
    Genotype = factor(as.character(Genotype),
                      levels = rev(c("trua1a", "trua1b", "trua2",
                                     "trua3", "trua4", "trua5", "trua6"))),
    Region = factor(as.character(Region), levels = c("5'UTR", "CDS", "3'UTR"))
  )

# methylKit can call zero dependent sites for some genotype × tissue panels.
# A 0/0 proportion is undefined, so drop those panels rather than aborting.
empty_panels <- summary_df |>
  distinct(Tissue, Genotype, Total) |>
  filter(Total == 0)
if (nrow(empty_panels) > 0) {
  message("[plot.R] Dropping ", nrow(empty_panels),
          " genotype × tissue panel(s) with zero dependent sites.")
  summary_df <- summary_df |> filter(Total > 0)
}
if (nrow(summary_df) == 0) {
  stop("No genotype/tissue has any annotated dependent sites; nothing to plot.")
}

write.table(
  summary_df,
  file.path(script_dir, "dependent_region_distribution_proportions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

fill_values <- PAL_REGION

file_suffix <- function(tissue_name) {
  gsub("[^A-Za-z0-9]+", "_", tolower(as.character(tissue_name)))
}

plot_region_distribution <- function(summary_data) {
  ggplot(summary_data, aes(x = Proportion, y = Genotype, fill = Region)) +
    geom_col(
      width = 0.72,
      color = "black",
      linewidth = base_linewidth,
      position = position_stack(reverse = TRUE)
    ) +
    scale_fill_manual(
      values = fill_values,
      breaks = c("5'UTR", "CDS", "3'UTR"),
      drop = TRUE,
      name = NULL
    ) +
    scale_x_continuous(
      labels = scales::percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0)),
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25)
    ) +
    scale_y_discrete(name = NULL, drop = TRUE) +
    labs(x = NULL, y = NULL) +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = "right",
      legend.direction = "vertical",
      legend.justification = "center",
      legend.key.spacing.y = unit(12, "points"),
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title = element_text(size = 22),
      legend.text = element_text(size = 18),
      plot.margin = margin(16, 20, 14, 18)
    )
}

plot_combined_region_distribution <- function(summary_data) {
  plot_region_distribution(summary_data) +
    facet_grid(Tissue ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(Tissue = c(Root = "Roots", Leaf = "Leaves", Flower = "Flowers"))) +
    theme(
      strip.background = element_rect(color = "black", fill = "white",
                                      linewidth = base_linewidth),
      strip.text.y = element_text(size = 20, angle = 0),
      axis.text.y = element_text(size = 18)
    )
}

combined_plot <- plot_combined_region_distribution(summary_df)
save_plot(
  combined_plot,
  file.path(script_dir, "dependent_region_distribution_barplot"),
  width = 11.5,
  height = 14.2
)

for (tissue_name in levels(droplevels(summary_df$Tissue))) {
  tissue_df <- summary_df |> filter(Tissue == tissue_name) |> select(-Tissue)

  write.table(
    tissue_df,
    file.path(script_dir, paste0("dependent_region_distribution_proportions_",
                                 file_suffix(tissue_name), ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  p <- plot_region_distribution(tissue_df)
  save_plot(
    p,
    file.path(script_dir, paste0("dependent_region_distribution_barplot_",
                                 file_suffix(tissue_name))),
    width = 11.5,
    height = 7.1
  )
}

message("[plot.R] Done: dependent region distribution barplots.")
print(summary_df)
