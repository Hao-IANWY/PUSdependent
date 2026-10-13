#!/usr/bin/env Rscript
# plot.R — WT modification-ratio raincloud plot with dependent sites highlighted.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(tibble)
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

input_file <- file.path(script_dir, "dependent_site_ratio_distribution_sites.tsv")
if (!file.exists(input_file)) stop("Input table not found: ", input_file, "\nRun data_prep.sh first.")

plot_df <- read.table(
  input_file,
  header = TRUE,
  sep = "\t",
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
) |>
  mutate(
    Tissue = factor(Tissue, levels = c("Root", "leaf", "Flower"),
                    labels = c("Root", "Leaf", "Flower")),
    WT_Ratio = as.numeric(WT_Ratio),
    Site_Status = factor(Site_Status, levels = c("wt_only", "dependent"),
                         labels = c("WT sites", "Dependent sites")),
    Tissue_Index = as.numeric(Tissue)
  ) |>
  filter(is.finite(WT_Ratio), WT_Ratio >= 0, WT_Ratio <= 1)

if (nrow(plot_df) == 0) stop("No valid WT sites found in: ", input_file)

summary_tbl <- plot_df |>
  group_by(Tissue) |>
  summarise(
    WT_Sites = n(),
    Dependent_Sites = sum(Site_Status == "Dependent sites"),
    WT_Sites_Above_20pct = sum(WT_Ratio >= 0.20),
    Dependent_Sites_Above_20pct = sum(WT_Ratio >= 0.20 & Site_Status == "Dependent sites"),
    WT_Median_Ratio = median(WT_Ratio),
    Dependent_Median_Ratio = median(WT_Ratio[Site_Status == "Dependent sites"]),
    WT_Fraction_Above_20pct = mean(WT_Ratio >= 0.20),
    Dependent_Fraction_Above_20pct = mean(WT_Ratio[Site_Status == "Dependent sites"] >= 0.20),
    .groups = "drop"
  )

write.table(
  summary_tbl,
  file.path(script_dir, "dependent_site_ratio_distribution_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

build_density_cloud <- function(site_tbl, max_width = 0.28, y_offset = 0.16) {
  tissue_levels <- levels(site_tbl$Tissue)
  lapply(seq_along(tissue_levels), function(i) {
    tissue_name <- tissue_levels[i]
    x <- site_tbl$WT_Ratio[site_tbl$Tissue == tissue_name]
    x <- x[is.finite(x)]
    if (length(unique(x)) < 2) return(tibble())
    den <- stats::density(x, from = min(x), to = max(x), n = 512, adjust = 1)
    width <- den$y / max(den$y, na.rm = TRUE) * max_width
    tibble(
      Tissue = factor(tissue_name, levels = tissue_levels),
      x = c(den$x, rev(den$x)),
      y = c(rep(i + y_offset, length(den$x)), rev(i + y_offset + width)),
      Cloud_ID = tissue_name
    )
  }) |> bind_rows()
}

cloud_tbl <- build_density_cloud(plot_df)
n_label_tbl <- summary_tbl |>
  mutate(
    Tissue_Index = as.numeric(Tissue),
    label = paste0(
      "n = ", scales::comma(WT_Sites), "; Dep = ", scales::comma(Dependent_Sites), "\n",
      ">=20%: ", scales::comma(Dependent_Sites_Above_20pct), "/", scales::comma(WT_Sites_Above_20pct), " Dep"
    ),
    x = 1.03
  )

tissue_colors <- setNames(PAL_TISSUE, c("Root", "Leaf", "Flower"))

set.seed(20260426)
p <- ggplot(plot_df, aes(x = WT_Ratio, y = Tissue_Index)) +
  geom_vline(xintercept = 0.20, linetype = "dashed", linewidth = base_linewidth,
             color = "grey35") +
  geom_polygon(
    data = cloud_tbl,
    aes(x = x, y = y, group = Cloud_ID, fill = Tissue),
    inherit.aes = FALSE,
    color = "black",
    linewidth = base_linewidth,
    alpha = 0.80
  ) +
  geom_boxplot(
    aes(y = Tissue_Index + 0.08, group = Tissue),
    width = 0.07,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = base_linewidth,
    orientation = "y"
  ) +
  geom_point(
    aes(y = Tissue_Index - 0.16, fill = "non-dependent sites"),
    position = position_jitter(width = 0, height = 0.06, seed = 20260426),
    shape = 21,
    size = 1.8,
    stroke = 0.25,
    color = "black",
    alpha = 0.45
  ) +
  geom_point(
    data = filter(plot_df, Site_Status == "Dependent sites"),
    aes(y = Tissue_Index - 0.16, fill = "Dependent sites"),
    position = position_jitter(width = 0, height = 0.06, seed = 20260426),
    shape = 21,
    size = 2.4,
    stroke = 0.35,
    color = "black",
    alpha = 0.95
  ) +
  geom_text(
    data = n_label_tbl,
    aes(x = x, y = Tissue_Index, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    fontface = "bold",
    size = 5.4,
    hjust = 0
  ) +
  scale_fill_manual(
    values = c(tissue_colors, setNames(PAL_STATUS, c("non-dependent sites", "Dependent sites"))),
    breaks = c("non-dependent sites", "Dependent sites"),
    name = NULL
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0)),
    name = "WT modification ratio"
  ) +
  scale_y_continuous(
    breaks = seq_along(levels(plot_df$Tissue)),
    labels = levels(plot_df$Tissue),
    expand = expansion(mult = c(0.10, 0.12)),
    name = NULL
  ) +
  coord_cartesian(xlim = c(0, 1.17), clip = "off") +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.text = element_text(size = 14),
    legend.key.size = unit(10, "points"),
    legend.key.spacing.y = unit(6, "points"),
    legend.margin = margin(4, 4, 4, 4),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20),
    axis.title.x = element_text(size = 22),
    plot.margin = margin(18, 92, 18, 18)
  )

save_plot(
  p,
  file.path(script_dir, "dependent_site_ratio_distribution_raincloud"),
  width = 10.4,
  height = 7
)
message("[plot.R] Done: dependent_site_ratio_distribution_raincloud.pdf/.png/.svg")
print(summary_tbl)
