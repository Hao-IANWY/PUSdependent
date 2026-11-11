#!/usr/bin/env Rscript
# plot.R — Mean RNAfold paired-probability profiles for four TRUA1A-dependent site groups
# Input: trua1a_groups_structure_unpaired_profile.tsv and trua1a_groups_structure_metrics.tsv
# Output: PDF/PNG/SVG profiles for the ±100-nt window, the Carlile-style
#         3′ (0–15 nt) region in downstream_15 and the extended −20 to +60 nt
#         window in window_-20_to_60 at the default 21 °C, plus a structure summary.
suppressPackageStartupMessages(library(ggplot2))

# Cairo replaces UTF-8 Greek letters and primes with ".." when R starts in the
# unavailable C.UTF-8 locale on this host. Use an installed UTF-8 locale before
# creating any text grobs so the title and axis labels keep Ψ and 3′ intact.
if (!nzchar(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))) {
  stop("Could not enable an installed UTF-8 locale (en_US.UTF-8).")
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) > 0) dirname(normalizePath(sub("^--file=", "", script_arg[1]))) else normalizePath(getwd())
shared_dir <- file.path(script_dir, "..", "..", "..", "..", "shared")
if (!file.exists(file.path(shared_dir, "plot_style.R"))) {
  stop("shared plot_style.R not found: ", file.path(shared_dir, "plot_style.R"))
}
source(file.path(shared_dir, "plot_style.R"))
source(file.path(shared_dir, "plot_tokens.R"))

read_tsv <- function(name) {
  f <- file.path(script_dir, name)
  if (!file.exists(f)) stop("Input table not found: ", f, "\nRun data_prep.sh first.")
  read.table(f, header = TRUE, sep = "\t", quote = "", comment.char = "",
             stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
}

profile <- read_tsv("trua1a_groups_structure_unpaired_profile.tsv")
metrics <- read_tsv("trua1a_groups_structure_metrics.tsv")
if (nrow(profile) == 0 || nrow(metrics) == 0) stop("Empty RNAfold input table.")

# Plant RNA structures are analyzed at the project-wide default temperature (21 °C).
default_temp_C <- 21
profile$temp_C <- as.numeric(profile$temp_C)
metrics$temp_C <- as.numeric(metrics$temp_C)
profile <- profile[profile$temp_C == default_temp_C, , drop = FALSE]
metrics <- metrics[metrics$temp_C == default_temp_C, , drop = FALSE]
if (nrow(profile) == 0 || nrow(metrics) == 0) {
  stop("No RNAfold rows found at the default temperature (21 °C).")
}

set_levels <- c("unique", "Root", "Leaf", "Flower")
set_labels <- c("Union (unique sites)", "Root-dependent", "Leaf-dependent", "Flower-dependent")
pal_group <- setNames(c("#6A7482CC", PAL_TISSUE_NPG[c("Root", "Leaf", "Flower")]), set_labels)
names(pal_group) <- set_labels

window_levels <- c(25, 50, 100)
window_labels <- paste0("±", window_levels, " nt")

profile$set <- factor(profile$set, levels = set_levels, labels = set_labels)
metrics$set <- factor(metrics$set, levels = set_levels, labels = set_labels)
profile$window <- factor(profile$window, levels = window_levels, labels = window_labels)
metrics$window <- factor(metrics$window, levels = window_levels, labels = window_labels)
profile$rel_pos <- as.numeric(profile$rel_pos)
profile$unpaired_prob <- as.numeric(profile$unpaired_prob)
# The folding output stores unpaired probability. Convert it to paired probability for all plots.
profile$paired_prob <- 1 - profile$unpaired_prob
metrics$center_unpaired_prob <- as.numeric(metrics$center_unpaired_prob)
metrics$mfe_per_nt <- as.numeric(metrics$mfe_per_nt)

sem <- function(x) if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else 0
prof_summary <- aggregate(paired_prob ~ set + window + temp_C + rel_pos,
                          data = profile, FUN = mean)
prof_sem <- aggregate(paired_prob ~ set + window + temp_C + rel_pos,
                      data = profile, FUN = sem)
names(prof_summary)[ncol(prof_summary)] <- "mean"
names(prof_sem)[ncol(prof_sem)] <- "sem"
prof_summary <- merge(prof_summary, prof_sem,
                      by = c("set", "window", "temp_C", "rel_pos"), sort = FALSE)

# Site counts are determined once per set/window (not once per temperature).
site_keys <- unique(metrics[, c("set", "window", "transcript", "pos", "strand")])
n_sites <- aggregate(cbind(n_sites = transcript) ~ set + window,
                     data = transform(site_keys, transcript = 1), FUN = sum)
set_n <- tapply(n_sites$n_sites, n_sites$set, max)
set_n <- set_n[set_labels]
legend_labels <- sprintf("%s (n=%d)", set_labels, as.integer(set_n))
names(legend_labels) <- set_labels
prof_summary$set <- factor(as.character(prof_summary$set), levels = set_labels,
                           labels = unname(legend_labels))
pal_plot <- setNames(unname(pal_group), unname(legend_labels))

init_theme(base_size = 16)

base_profile_plot <- function(df, title = NULL) {
  p <- ggplot(df, aes(x = rel_pos, y = mean, color = set, fill = set)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.14, color = NA) +
    geom_line(linewidth = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    scale_color_manual(values = pal_plot, name = NULL) +
    scale_fill_manual(values = pal_plot, name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Position relative to Ψ site (nt)", y = "Mean paired probability", title = title) +
    theme(legend.position = "top", axis.text.x = element_text(size = 9))
  if (!is.null(title)) p <- p + theme(plot.title = element_text(size = 12))
  p
}

# Keep the focused-region deliverables in their own subfolders.
region_plots <- list(
  list(dir = "downstream_15", suffix = "3prime", lo = 0, hi = 15,
       breaks = seq(0, 15, 5), width = 7.2, points = TRUE,
       title = "Paired 3′ of Ψ"),
  list(dir = "window_-20_to_60", suffix = "-20_to_60", lo = -20, hi = 60,
       breaks = seq(-20, 60, 10), width = 9.5, points = FALSE,
       title = "Paired around Ψ (−20 to +60 nt)")
)

{
  # Temperature is fixed at 21 °C and is intentionally not shown in figure filenames or labels.
  sub <- prof_summary

  group_guide <- guide_legend(nrow = 2, byrow = TRUE)
  p_100 <- base_profile_plot(sub[sub$window == "±100 nt", ]) +
    scale_x_continuous(breaks = seq(-100, 100, 25)) +
    guides(colour = group_guide, fill = group_guide)
  save_plot(p_100, file.path(script_dir, "trua1a_groups_paired_profile_100nt"),
            width = 8, height = 4.8)

  for (reg in region_plots) {
    reg_dir <- file.path(script_dir, reg$dir)
    dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
    reg_data <- sub[sub$window == "±100 nt" & sub$rel_pos >= reg$lo & sub$rel_pos <= reg$hi, ]
    p_reg <- base_profile_plot(reg_data, title = reg$title) +
      scale_x_continuous(breaks = reg$breaks) +
      guides(colour = group_guide, fill = group_guide) +
      theme(plot.title = element_text(size = 13))
    if (reg$points) p_reg <- p_reg + geom_point(size = 1.7)
    save_plot(p_reg, file.path(reg_dir, paste0("trua1a_groups_", reg$suffix, "_paired")),
              width = reg$width, height = 5.0)
  }
}

# Compact center-site and folding summary.
center_summary <- do.call(rbind, lapply(
  split(metrics, list(metrics$set, metrics$window, metrics$temp_C), drop = TRUE),
  function(df) {
    data.frame(
      set = as.character(df$set[1]),
      window = as.character(df$window[1]),
      temp_C = df$temp_C[1],
      n_sites = length(unique(paste(df$transcript, df$pos, df$strand))),
      mean_center_unpaired_prob = round(mean(df$center_unpaired_prob), 4),
      median_mfe_per_nt = round(stats::median(df$mfe_per_nt), 4),
      frac_center_paired_MFE = round(mean(df$center_mfe_state == "paired"), 4)
    )
  }
))
write.table(center_summary, file.path(script_dir, "trua1a_groups_structure_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[plot] Wrote: ", file.path(script_dir, "trua1a_groups_structure_summary.tsv"))
message("[plot] Done: four-group RNAfold paired-probability figures")
