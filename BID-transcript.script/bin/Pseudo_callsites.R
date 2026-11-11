# Calls candidate pseudouridine (Ψ) sites from cpup-format deletion-signal
# matrices using a replicate-aware one-sided Fisher's exact test.
#
# Inputs:  per-replicate TSVs with columns chr, pos, strand, context, ratio, depth, gap
# Outputs: <output> full results table; <output>.positive.tsv passing sites
#          (chr, pos, strand, mean Treat ratio/depth over passing replicates)
#
# Calling criteria:
#   1. Mean Mock and Treat depth >= --min-avg-depth (default: 10).
#   2. Significant one-sided Fisher's test on counts pooled across
#      biological replicates (BH-adjusted p < --max-padj, default: 0.01).
#   3. Per-replicate filters:
#      - Treat gap count >= --min-gap (default: 5)
#      - Treat gap ratio >= --min-treat-ratio (default: 0.05)
#      - Treat gap ratio >= ratio_fold * Mock gap ratio (default: 2x)
#
# A site is positive when at least --min-pass-reps replicates pass.

library(dplyr)
library(parallel)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input files (e.g., Mock:Mock,Treat:Treat; Mock must precede Treat)", dest = "input"),
  make_option(c("-n", "--name"), type = "character", help = "Sample names (e.g., Mock:Mock,Treat:Treat; Mock must precede Treat)", dest = "name"),
  make_option(c("-k", "--keepall"), action = "store_true", default = FALSE,
              help = "Keep all site information (default: FALSE)", dest = "keepall"),
  make_option(c("-r", "--min-treat-ratio"), type = "double", default = 0.05,
              help = "Minimum Treat ratio threshold per replicate (default: 0.05)", dest = "min_treat_ratio"),
  make_option(c("--min-avg-depth"), type = "double", default = 10,
              help = "Minimum Mock/Treat mean depth (default: 10)", dest = "min_avg_depth"),
  make_option(c("--max-padj"), type = "double", default = 0.01,
              help = "BH-adjusted p-value cutoff (default: 0.01)", dest = "max_padj"),
  make_option(c("--min-gap"), type = "double", default = 5,
              help = "Minimum gap count threshold per replicate (default: 5)", dest = "min_gap"),
  make_option(c("--ratio-fold"), type = "double", default = 2,
              help = "Required fold increase: Treat ratio >= fold * Mock ratio (default: 2)", dest = "ratio_fold"),
  make_option(c("--min-pass-reps"), type = "integer", default = 2,
              help = "Minimum replicates required for site to be considered positive (default: 2)", dest = "min_pass_reps"),
  make_option(c("-c", "--read-cores"), type = "integer", default = max(1, detectCores() - 1),
              help = "Cores for file reading (default: detectCores() - 1)", dest = "read_cores"),
  make_option(c("-o", "--output"), type = "character", help = "Output file path", dest = "output")
)

# Parse command-line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# Validate arguments
if (is.null(opt$input) || is.null(opt$name) || is.null(opt$output)) {
  stop("Please provide all required arguments: -i for input, -n for name, -o for output")
}
if (is.na(opt$min_treat_ratio) || opt$min_treat_ratio < 0) {
  stop("--min-treat-ratio must be a non-negative number")
}
if (is.na(opt$min_avg_depth) || opt$min_avg_depth <= 0) {
  stop("--min-avg-depth must be a positive number")
}
if (is.na(opt$max_padj) || opt$max_padj <= 0 || opt$max_padj > 1) {
  stop("--max-padj must be in (0, 1]")
}
if (is.na(opt$min_gap) || opt$min_gap < 0) {
  stop("--min-gap must be a non-negative number")
}
if (is.na(opt$ratio_fold) || opt$ratio_fold <= 0) {
  stop("--ratio-fold must be a positive number")
}
if (is.na(opt$min_pass_reps) || opt$min_pass_reps < 1) {
  stop("--min-pass-reps must be a positive integer")
}
if (is.na(opt$read_cores) || opt$read_cores < 1) {
  stop("--read-cores must be a positive integer")
}

Input <- opt$input
Input_names <- opt$name
output <- opt$output
keepall <- opt$keepall
min_treat_ratio <- opt$min_treat_ratio
min_avg_depth <- opt$min_avg_depth
max_padj <- opt$max_padj
min_gap <- opt$min_gap
ratio_fold <- opt$ratio_fold
min_pass_reps <- opt$min_pass_reps
read_cores <- opt$read_cores
positive_output <- if (grepl("\\.Fisher\\.tsv$", output)) {
  sub("\\.Fisher\\.tsv$", ".positive.tsv", output)
} else if (grepl("\\.tsv$", output)) {
  sub("\\.tsv$", ".positive.tsv", output)
} else {
  paste0(output, ".positive.tsv")
}

message("========== Parameter Settings ==========")
message("Input files                 : ", Input)
message("Sample names                : ", Input_names)
message("Output                      : ", output)
message("Positive output             : ", positive_output)
message("Keep all sites              : ", keepall)
message("Min Treat ratio             : ", min_treat_ratio)
message("Min avg depth               : ", min_avg_depth)
message("Max p.adj                   : ", max_padj)
message("Min gap                     : ", min_gap)
message("Ratio fold                  : ", ratio_fold)
message("Min pass reps               : ", min_pass_reps)
message("Read cores                  : ", read_cores)
message("========================================")

# Mock entries must precede Treat entries in --input and --name
file_group <- strsplit(Input, ",")[[1]]
files <- unlist(lapply(file_group, function(x) strsplit(x, ":")[[1]]))

# Extract sample names and the Mock/Treat name lists
names_group <- strsplit(Input_names, ",")[[1]]
sample_names <- unlist(lapply(names_group, function(x) strsplit(x, ":")[[1]]))
mock_sample_names <- strsplit(names_group[1], ":")[[1]]
treat_sample_names <- strsplit(names_group[2], ":")[[1]]

# File reading utility
read_file <- function(file, sample_name) {
  # Check if file exists
  if (!file.exists(file)) {
    stop(sprintf("Error: File not found: %s (Sample: %s)", file, sample_name))
  }
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data <- data[, c("chr", "pos", "strand", "context", "ratio", "depth", "gap")]
  # Enforce numeric types
  numeric_cols <- c("ratio", "depth", "gap")
  data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)
  # Rename columns to sample-specific names
  colnames(data)[5:7] <- paste0(sample_name, c("_ratio", "_depth", "_gap"))
  return(data)
}
# Read all files
data_list <- mclapply(seq_along(files), function(i) 
  read_file(files[i], sample_names[i]), mc.cores = read_cores)
# Merge based on site columns
merged_data <- Reduce(function(x, y) merge(x, y, by = c("chr", "pos", "strand", "context"), all = TRUE), data_list)
# Fill missing values with 0
numeric_cols <- sapply(merged_data, is.numeric)
merged_data[numeric_cols][is.na(merged_data[numeric_cols])] <- 0

# Calculate group-wise averages
Mock_depth <- paste0(strsplit(names_group[1], ":")[[1]], "_depth")
Treat_depth <- paste0(strsplit(names_group[2], ":")[[1]], "_depth")
Mock_ratio <- paste0(strsplit(names_group[1], ":")[[1]], "_ratio")
Treat_ratio <- paste0(strsplit(names_group[2], ":")[[1]], "_ratio")
Mock_gap <- paste0(strsplit(names_group[1], ":")[[1]], "_gap")
Treat_gap <- paste0(strsplit(names_group[2], ":")[[1]], "_gap")
# Compute average depths
merged_data$Mock_avg_depth <- rowMeans(as.data.frame(merged_data[, Mock_depth]))
merged_data$Treat_avg_depth <- rowMeans(as.data.frame(merged_data[, Treat_depth]))
# Compute average ratios
merged_data$Mock_avg_ratio <- rowMeans(as.data.frame(merged_data[, Mock_ratio]))
merged_data$Treat_avg_ratio <- rowMeans(as.data.frame(merged_data[, Treat_ratio]))
# Compute average gaps
merged_data$Treat_avg_gap <- rowMeans(as.data.frame(merged_data[, Treat_gap]))

# Filter if !keepall
if (!keepall) {
	# Remove sites with depth below min_avg_depth
	merged_data <- merged_data[merged_data$Mock_avg_depth >= min_avg_depth & merged_data$Treat_avg_depth >= min_avg_depth, ]
}

# One-sided Fisher's Exact Test
di <- rowSums(merged_data[Mock_depth])
gi <- rowSums(merged_data[Mock_gap])
dt <- rowSums(merged_data[Treat_depth])
gt <- rowSums(merged_data[Treat_gap])
get_p <- function(r) {
  fisher.test(as.table(matrix(as.numeric(r[c(1:4)]), ncol = 2)), alt = "greater")$p.value
}

merged_data$p.value <- apply(cbind(di - gi, gi, dt - gt, gt), 1, get_p)
merged_data$p.adj <- p.adjust(merged_data$p.value, method = "BH")
merged_data <- merged_data[merged_data$p.adj < max_padj, ]

# Perform per-replicate pass assessment
# Pairwise logic for same counts; mean-based logic for different counts.
if (length(mock_sample_names) == length(treat_sample_names)) {
  pass_treat_ratio <- Treat_ratio
  pass_treat_depth <- Treat_depth

  for (i in seq_along(mock_sample_names)) {
    mock_name <- mock_sample_names[i]
    treat_name <- treat_sample_names[i]
    pass_col <- paste0("Rep", i, "_pass")
    pvalue_col <- paste0("Rep", i, "_pvalue")
    padj_col <- paste0("Rep", i, "_p.adj")

    # Perform Fisher's test
    mock_depth <- merged_data[paste0(mock_name, "_depth")]
    mock_gap <- merged_data[paste0(mock_name, "_gap")]
    treat_depth <- merged_data[paste0(treat_name, "_depth")]
    treat_gap <- merged_data[paste0(treat_name, "_gap")]
    merged_data[[pvalue_col]] <- apply(cbind(mock_depth - mock_gap, mock_gap,
                                              treat_depth - treat_gap, treat_gap), 1, get_p)
    merged_data[[padj_col]] <- p.adjust(merged_data[[pvalue_col]], method = "BH")

    merged_data[[pass_col]] <- ifelse(
      as.numeric(merged_data[[paste0(treat_name, "_gap")]]) >= min_gap &
      as.numeric(merged_data[[paste0(treat_name, "_ratio")]]) >= min_treat_ratio &
      as.numeric(merged_data[[paste0(treat_name, "_ratio")]]) >= ratio_fold * (as.numeric(merged_data[[paste0(mock_name, "_ratio")]]) + 0.0001),
      1,
      0
    )
  }
} else {
  message("Mock/Treat replicate counts differ; using Mock group mean for per-Treat pass calls.")
  merged_data$Mock_mean_depth <- rowMeans(as.data.frame(merged_data[, Mock_depth, drop = FALSE]))
  merged_data$Mock_mean_ratio <- rowMeans(as.data.frame(merged_data[, Mock_ratio, drop = FALSE]))
  merged_data$Mock_mean_gap <- rowMeans(as.data.frame(merged_data[, Mock_gap, drop = FALSE]))

  pass_treat_ratio <- Treat_ratio
  pass_treat_depth <- Treat_depth

  for (i in seq_along(treat_sample_names)) {
    treat_name <- treat_sample_names[i]
    pass_col <- paste0("TreatRep", i, "_pass")
    pvalue_col <- paste0("TreatRep", i, "_pvalue")
    padj_col <- paste0("TreatRep", i, "_p.adj")

    mock_depth <- merged_data["Mock_mean_depth"]
    mock_gap <- merged_data["Mock_mean_gap"]
    treat_depth <- merged_data[paste0(treat_name, "_depth")]
    treat_gap <- merged_data[paste0(treat_name, "_gap")]
    merged_data[[pvalue_col]] <- apply(cbind(mock_depth - mock_gap, mock_gap,
                                              treat_depth - treat_gap, treat_gap), 1, get_p)
    merged_data[[padj_col]] <- p.adjust(merged_data[[pvalue_col]], method = "BH")

    merged_data[[pass_col]] <- ifelse(
      as.numeric(merged_data[[paste0(treat_name, "_gap")]]) >= min_gap &
      as.numeric(merged_data[[paste0(treat_name, "_ratio")]]) >= min_treat_ratio &
      as.numeric(merged_data[[paste0(treat_name, "_ratio")]]) >= ratio_fold * (as.numeric(merged_data$Mock_mean_ratio) + 0.0001),
      1,
      0
    )
  }
}

# Generate total pass columns
pass_columns <- grep("_pass$", colnames(merged_data), value = TRUE)
# Final status: pass if >= min_pass_reps replicates pass
merged_data$Pass <- as.integer(rowSums(merged_data[pass_columns]) >= min_pass_reps)

# Export results
positive_data <- merged_data[merged_data$Pass == 1, c("chr", "pos", "strand"), drop = FALSE]

if (nrow(positive_data) > 0) {
  pass_matrix <- as.matrix(merged_data[merged_data$Pass == 1, pass_columns, drop = FALSE])
  treat_ratio_matrix <- as.matrix(merged_data[merged_data$Pass == 1, pass_treat_ratio, drop = FALSE])
  treat_depth_matrix <- as.matrix(merged_data[merged_data$Pass == 1, pass_treat_depth, drop = FALSE])
  pass_counts <- rowSums(pass_matrix)

  positive_data$ratio <- rowSums(treat_ratio_matrix * pass_matrix) / pass_counts
  positive_data$depth <- rowSums(treat_depth_matrix * pass_matrix) / pass_counts
} else {
  positive_data$ratio <- numeric(0)
  positive_data$depth <- numeric(0)
}

write.table(merged_data, output, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(positive_data, positive_output, sep = "\t", row.names = FALSE, quote = FALSE)
