#!/usr/bin/env Rscript
# plot_union_paired_profile.R — Mean paired-probability profile of the Union
# TRUA1A-dependent Ψ sites across the ±100-nt window and the downstream 0–15 nt
# region at the default 21 °C. Site count is read from the profile table.
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

  # Focused Carlile-style downstream region, using the same single-group styling
  # as the ±100-nt Union profile above.
  downstream_dir <- file.path(script_dir, "downstream_15")
  dir.create(downstream_dir, recursive = TRUE, showWarnings = FALSE)
  p_downstream <- ggplot(prof[prof$rel_pos >= 0 & prof$rel_pos <= 15, ],
                         aes(x = rel_pos, y = mean)) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), alpha = 0.18, fill = col_red) +
    geom_line(linewidth = 0.75, color = col_red) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 15)) +
    scale_x_continuous(breaks = seq(0, 15, 5)) +
    labs(x = "Position relative to Ψ site (nt)", y = "Mean paired probability",
         title = sprintf("Union TRUA1A-dependent Ψ sites (n=%d)", n_sites)) +
    theme(plot.title = element_text(size = 12))
  save_plot(p_downstream, file.path(downstream_dir, "trua1a_union_3prime_paired"),
            width = 7.2, height = 5.0)
}

message(sprintf("[plot] Done: Union %d-site paired-probability profiles (21 °C; not annotated in figures)", n_sites))
