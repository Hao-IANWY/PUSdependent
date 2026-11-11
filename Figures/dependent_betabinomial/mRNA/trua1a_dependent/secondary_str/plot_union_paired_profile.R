#!/usr/bin/env Rscript
# plot_union_paired_profile.R — Mean paired-probability profile of the Union
# TRUA1A-dependent Ψ sites across the ±100-nt window, the downstream 0–15 nt
# region (downstream_15) and the extended −20 to +60 nt window
# (window_-20_to_60) at the default 21 °C. Site count is read from the profile table.
# Input: trua1a_groups_structure_unpaired_profile.tsv (run data_prep.sh first)
suppressPackageStartupMessages(library(ggplot2))

# Keep Unicode labels intact in Cairo outputs when the shell starts in C.UTF-8.
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

profile <- read.table(file.path(script_dir, "trua1a_groups_structure_unpaired_profile.tsv"),
                      header = TRUE, sep = "\t", quote = "", comment.char = "",
                      stringsAsFactors = FALSE, check.names = FALSE, na.strings = "")
sub <- profile[profile$set == "unique" & profile$window == 100 &
               as.numeric(profile$temp_C) == 21, ]
if (nrow(sub) == 0) stop("No Union ±100-nt rows found at 21 °C.")
n_sites <- nrow(unique(sub[, c("transcript", "pos")]))
sub$unpaired_prob <- as.numeric(sub$unpaired_prob)
sub$rel_pos <- as.numeric(sub$rel_pos)
sub$paired_prob <- 1 - sub$unpaired_prob

sem <- function(x) if (length(x) > 1) stats::sd(x) / sqrt(length(x)) else 0
prof_mean <- aggregate(paired_prob ~ rel_pos, data = sub, FUN = mean)
prof_sem <- aggregate(paired_prob ~ rel_pos, data = sub, FUN = sem)
names(prof_mean)[2] <- "mean"
names(prof_sem)[2] <- "sem"
prof <- merge(prof_mean, prof_sem, by = "rel_pos", sort = FALSE)

init_theme(base_size = 16)

col_red <- "#E64B35"  # NPG red, project-wide "dependent" color

{
  p <- ggplot(prof, aes(x = rel_pos, y = mean)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, fill = col_red) +
    geom_line(linewidth = 0.75, color = col_red) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(breaks = seq(-100, 100, 25)) +
    labs(x = "Position relative to Ψ site (nt)", y = "Mean paired probability",
         title = sprintf("Union TRUA1A-dependent Ψ sites (n=%d)", n_sites)) +
    theme(plot.title = element_text(size = 12))
  save_plot(p, file.path(script_dir, "trua1a_union_paired_profile_100nt"),
            width = 8, height = 4.6)

  # Focused regions, using the same single-group styling as the ±100-nt Union
  # profile above: the Carlile-style 3′ window and the extended −20 to +60 nt window.
  region_plots <- list(
    list(dir = "downstream_15", suffix = "3prime", lo = 0, hi = 15,
         breaks = seq(0, 15, 5), width = 7.2),
    list(dir = "window_-20_to_60", suffix = "-20_to_60", lo = -20, hi = 60,
         breaks = seq(-20, 60, 10), width = 9.5)
  )
  for (reg in region_plots) {
    reg_dir <- file.path(script_dir, reg$dir)
    dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)
    p_reg <- ggplot(prof[prof$rel_pos >= reg$lo & prof$rel_pos <= reg$hi, ],
                    aes(x = rel_pos, y = mean)) +
      geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, fill = col_red) +
      geom_line(linewidth = 0.75, color = col_red) +
      geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
      coord_cartesian(ylim = c(0, 1), xlim = c(reg$lo, reg$hi)) +
      scale_x_continuous(breaks = reg$breaks) +
      labs(x = "Position relative to Ψ site (nt)", y = "Mean paired probability",
           title = sprintf("Union TRUA1A-dependent Ψ sites (n=%d)", n_sites)) +
      theme(plot.title = element_text(size = 12))
    save_plot(p_reg, file.path(reg_dir, paste0("trua1a_union_", reg$suffix, "_paired")),
              width = reg$width, height = 5.0)
  }
}

message(sprintf("[plot] Done: Union %d-site paired-probability profiles (21 °C; not annotated in figures)", n_sites))
