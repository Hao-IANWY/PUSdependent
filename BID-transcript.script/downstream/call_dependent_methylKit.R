#!/usr/bin/env Rscript
# ============================================================================
# call_dependent_methylKit.R - identify writer-dependent Ψ sites with methylKit
#
# Companion to call_dependent_betabinomial.R (glmmTMB beta-binomial version);
# this script instead tests the modification difference with methylKit,
# following the canonical workflow of Akalin et al. 2012 (Genome Biology
# 13:R87):
#
#   methRead()           read per-sample modification counts
#   unite()              merge samples (complete-case: a site must appear in
#                        every sample)
#   calculateDiffMeth()  differential-modification test. methylKit dispatches
#                        on the total sample count: with exactly 2 samples
#                        (1 per group) the test is forced to fast.fisher
#                        (Fisher's exact test); with biological replicates
#                        (3v3 in this project) it uses logReg, a binomial
#                        logistic regression weighted by coverage.
#                        overdispersion="MN" estimates the overdispersion
#                        factor phi; phi>1 triggers an F test
#                        (deviance/phi ~ F(1, n-2)), and phi<=1 degrades to
#                        chi-squared.
#
# Grouping convention: WT is treatment=1 and the mutant is treatment=0, so
#     meth.diff = (WT coverage-weighted mean modification) - (mutant
# coverage-weighted mean modification), in percentage points.
# Writer-dependent sites (modification lost in the mutant) are the
# meth.diff > 0 direction.
#
# Significance criteria (--significance):
#   q (default)  use the calculateDiffMeth qvalue (--methylkit-adjust, BH by
#                default) plus a positive direction (meth.diff > 0).
#   p            raw one-sided p (p_two/2, or 1 for the opposite direction),
#                retained for reproducing p-based results.
# The qvalue is a two-sided FDR. BH is the default adjustment; SLIM's
# per-line pi0 estimate is systematically underestimated on this project's
# 100-200 site candidate universe, giving liberal or even negative q values,
# so --methylkit-adjust SLIM is retained only for reproducing SLIM-based
# results. The one-sided p value is still written to the output p.value
# column for auditing and downstream reference.
#
# Output conventions follow call_dependent_betabinomial.R / call_dependent.py:
#   WT_Ratio  = --ratio-col (mean WT Treat modification fraction)
#   Mut_Ratio = simple mean of the mutant Treat replicate ratios
#   Delta_Ratio(WT-Mut) = WT_Ratio - Mut_Ratio (simple difference of means,
#                used for calling)
# methylKit's meth.diff is used only for direction and significance, never
# for delta. Thresholds: --alpha 0.05 and --min-delta 0.05 (fraction); a site
# is dependent only when both mutant lines pass.
#
# Depth filter: every mutant Mock and Treat replicate, and every WT replicate,
# must have depth >= --min-depth (default 10); a site failing any replicate is
# dropped (status="depth_filter"). The adjusted full-site grid writes depth=0
# and ratio=0 for uncovered sites (a division-guard convention, not a valid
# "modification lost" observation), so this filter keeps any zero- or
# low-coverage replicate out of the test and no zero-coverage special case is
# needed.
# ============================================================================
# Method overview
# ---------------------------------------------------------------------------
# Rationale. If a Ψ site is written by a specific PUS enzyme, mutating that
# enzyme should abolish the modification: the site's gap fraction in Treat
# libraries should fall back to the Mock background. The script therefore
# asks one question per site: is the mutant modification fraction
# significantly lower than WT? It answers separately for each of the two
# independent mutant lines of the genotype.
#
# Why a DNA-methylation toolkit applies. At each site and replicate, BID-seq
# yields exactly two numbers: how many reads cover the site and how many of
# them carry a deletion. This maps one-to-one onto a BS-seq cytosine record:
# coverage depth is sequencing coverage, and the gap count is the methylated-
# read count. Treating the gap fraction as the "methylation frequency" turns
# the test into a standard differential-methylation problem that methylKit
# solves directly.
#
# Per mutant line:
#
#   1. Candidate pool. Only sites called positive in WT of the same dataset
#      are tested; in the mutant we only ask whether WT-modified sites stay
#      modified.
#   2. Depth filter. Every WT Treat replicate and every mutant Mock and Treat
#      replicate must have depth >= 10 at the site; any failing replicate
#      drops the site. This prevents "not sequenced" from masquerading as
#      "modification lost" (uncovered sites carry depth = 0 in the adjusted
#      tables) and keeps only fully comparable sites. Mock libraries witness
#      coverage only; the compared signal is the Treat gap fraction.
#   3. Test. Per-site coverage and gap counts of WT and mutant replicates go
#      to methylKit (methRead -> unite -> calculateDiffMeth) for a group
#      comparison: coverage-weighted testing of whether the gap fraction
#      changes with genotype (WT versus mutant); the model also absorbs
#      replicate-level variation beyond binomial expectation
#      (overdispersion), so normal fluctuation is not mistaken for
#      significance. Hundreds to thousands of sites are tested at once, so
#      two-sided p values are BH-corrected into q values (two-sided FDR) to
#      control false positives across sites; calling uses that q value
#      directly. Because only the "mutant loses modification" direction is
#      of interest, the direction must also agree (meth.diff > 0, WT more
#      modified); the opposite direction never passes.
#   4. Effect size. Delta = WT_Ratio - Mut_Ratio, each term the arithmetic
#      mean of per-replicate Treat gap fractions. It is computed directly
#      from the input tables and independently of methylKit, so the reported
#      difference always matches the raw fractions.
#   5. Call. BH q < 0.05, positive direction (meth.diff > 0), and
#      Delta >= 0.05 make the site dependent in that line. The final call
#      keeps only sites dependent in both lines, so a line-specific
#      background cannot carry a site; the reported mutant ratio averages
#      the two lines.
#
# Outputs: per-line <allele>.dependent-sites.tsv, merged <gene>.dependent.tsv,
# and a summary JSON.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
})

# Load methylKit lazily; the main workflow attaches it when needed, so a
# missing package does not break --help.
if (!requireNamespace("methylKit", quietly = TRUE)) {
  stop("This script requires the 'methylKit' package. Install with:\n",
       "  conda install -c bioconda bioconductor-methylkit")
}

option_list <- list(
  make_option(c("--tsv"), type = "character",
              help = "WT candidate TSV with a header and coordinate + ratio columns."),
  make_option(c("--key-cols"), type = "character", default = "chr,pos,strand",
              help = "Comma-separated site key columns in adjusted TSV files. Default: chr,pos,strand."),
  make_option(c("--tsv-key-cols"), type = "character", default = NULL,
              help = "Comma-separated site key columns in --tsv. Default: use --key-cols."),
  make_option(c("--context-col"), type = "character", default = "context",
              help = "Optional context column name. Default: context."),
  make_option(c("--ratio-col"), type = "character", default = "Treat_avg_ratio",
              help = "WT ratio column in --tsv. Default: Treat_avg_ratio."),
  make_option(c("--wt-ratio-cols"), type = "character", default = NULL,
              help = "Comma-separated WT per-replicate Treat ratio columns. Default: auto-detect '^Treat-.*_ratio$'."),
  make_option(c("--wt-depth-cols"), type = "character", default = NULL,
              help = "Comma-separated WT per-replicate Treat depth columns. Default: auto-detect '^Treat-.*_depth$'."),
  make_option(c("--genotype1"), type = "character",
              help = "Name of the first mutant line, e.g. trua1a_1-Leaf."),
  make_option(c("--genotype1-treat-files"), type = "character",
              help = "Comma-separated Treat adjusted TSV files for genotype1."),
  make_option(c("--genotype1-mock-files"), type = "character",
              help = "Comma-separated Mock adjusted TSV files for genotype1."),
  make_option(c("--genotype2"), type = "character",
              help = "Name of the second mutant line, e.g. trua1a_2-Leaf."),
  make_option(c("--genotype2-treat-files"), type = "character",
              help = "Comma-separated Treat adjusted TSV files for genotype2."),
  make_option(c("--genotype2-mock-files"), type = "character",
              help = "Comma-separated Mock adjusted TSV files for genotype2."),
  make_option(c("--outdir"), type = "character",
              help = "Output directory for line-level, final-level, and summary files."),
  make_option(c("--output-prefix"), type = "character", default = NULL,
              help = "Prefix for final output files. Default: inferred from genotype1/2."),
  make_option(c("--adjusted-ratio-col"), type = "character", default = "ratio",
              help = "Ratio column in adjusted TSV files. Default: ratio."),
  make_option(c("--adjusted-depth-col"), type = "character", default = "depth",
              help = "Depth column in adjusted TSV files. Default: depth."),
  make_option(c("--min-depth"), type = "double", default = 10,
              help = "Minimum per-replicate depth required for every Mock/Treat (and WT) replicate. Default: 10."),
  make_option(c("--min-delta"), type = "double", default = 0.05,
              help = "Effect-size floor (WT_Ratio - Mut_Ratio), in fraction. Default: 0.05."),
  make_option(c("--alpha"), type = "double", default = 0.05,
              help = "Significance cutoff for --significance. Default: 0.05."),
  make_option(c("--significance"), type = "character", default = "q",
              help = "Statistic used for dependent calls: q = methylKit qvalue (default), p = raw one-sided p."),
  make_option(c("--methylkit-test"), type = "character", default = "F",
              help = "methylKit test statistic for the multi-replicate logistic-regression path: F (default; overdispersion-corrected F test, degrades to Chisq when phi<=1) or Chisq. Only used when total samples > 2; with exactly 2 samples methylKit forces fast.fisher regardless."),
  make_option(c("--methylkit-adjust"), type = "character", default = "BH",
              help = "Multiple-testing adjustment in calculateDiffMeth: BH (default) or SLIM (its per-run pi0 estimate is unreliable on small candidate sets)."),
  make_option(c("--methylkit-context"), type = "character", default = "CpH",
              help = "methylKit context label (metadata only for generic read). Default: CpH."),
  make_option(c("--assembly"), type = "character", default = "custom",
              help = "methylKit assembly label (metadata only). Default: custom."),
  make_option(c("--mc-cores"), type = "integer", default = 1,
              help = "Number of cores passed to methylKit's unite/calculateDiffMeth. Default: 1."),
  make_option(c("--allow-duplicate-sites"), action = "store_true", default = FALSE,
              help = "Keep first occurrence of duplicate keys instead of failing.")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("tsv", "genotype1", "genotype1-treat-files", "genotype1-mock-files",
              "genotype2", "genotype2-treat-files", "genotype2-mock-files", "outdir")
for (r in required) {
  if (is.null(opt[[r]]) || !nzchar(opt[[r]])) stop(sprintf("Missing required option: --%s", r))
}
if (!opt[["methylkit-test"]] %in% c("F", "Chisq")) {
  stop("--methylkit-test must be 'F' or 'Chisq'")
}
if (!opt[["methylkit-adjust"]] %in% c("BH", "SLIM")) {
  stop("--methylkit-adjust must be 'BH' or 'SLIM'")
}
if (!opt[["significance"]] %in% c("p", "q")) {
  stop("--significance must be 'p' or 'q'")
}
if (!(opt[["min-delta"]] >= 0 && opt[["min-delta"]] <= 1)) {
  stop("--min-delta must be a fraction in [0, 1]")
}

split_csv <- function(x) {
  items <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  items[nzchar(items)]
}

key_cols <- split_csv(opt[["key-cols"]])
tsv_key_cols <- if (is.null(opt[["tsv-key-cols"]]) || !nzchar(opt[["tsv-key-cols"]])) key_cols else split_csv(opt[["tsv-key-cols"]])
if (length(tsv_key_cols) != length(key_cols)) {
  stop("--tsv-key-cols and --key-cols must define the same number of columns")
}

make_key <- function(dt, cols) {
  do.call(paste, c(as.list(dt)[cols], list(sep = "|")))
}

# Read adjusted TSV, keeping only (ratio, depth) for target sites.
read_adjusted <- function(path, target_keys, allow_dup) {
  # The adjusted file can contain tens of millions of rows; pre-filter with grep -F
  # using (chr, pos, strand) to avoid full ingestion, passing only target site rows to fread.
  # (The "\t" in the key anchors pos/strand boundaries; chr is a unique transcript ID.)
  key_file <- tempfile(fileext = ".keys")
  writeLines(gsub("|", "\t", target_keys, fixed = TRUE), key_file)
  hdr <- strsplit(readLines(path, n = 1), "\t")[[1]]
  cmd <- sprintf("LC_ALL=C grep -F -f %s %s", shQuote(key_file), shQuote(path))
  dt <- fread(cmd = cmd, header = FALSE, col.names = hdr, showProgress = FALSE)
  setnames(dt, opt[["adjusted-ratio-col"]], "ratio")
  setnames(dt, opt[["adjusted-depth-col"]], "depth")
  dt[, site_id := make_key(.SD, key_cols)]
  dt <- dt[site_id %in% target_keys]
  dt[, ratio := as.numeric(ratio)]
  dt[, depth := as.numeric(depth)]
  if (!allow_dup && anyDuplicated(dt$site_id)) {
    stop(sprintf("%s: duplicate site key", path))
  }
  dt[, .(site_id, ratio, depth)]
}

# Write WT/Mutant samples to methylKit generic format (chr start end strand coverage freqC).
# freqC uses fraction (0-1); methylKit internally sets numCs = round(coverage * freqC).
write_methyl_sample <- function(path, key_dt, coverage, freqC) {
  out <- data.table(
    chr     = key_dt[[key_cols[1]]],
    start   = key_dt[[key_cols[2]]],
    end     = key_dt[[key_cols[2]]],
    strand  = key_dt[[key_cols[3]]],
    coverage = as.numeric(coverage),
    freqC    = as.numeric(freqC)
  )
  out[is.na(freqC), freqC := 0]
  fwrite(out, path, sep = "\t", quote = FALSE, na = "0")
}

# Column mapping for generic methylKit reading
pipeline_spec <- list(
  fraction = TRUE,   # freqC is a fraction (0-1), not 0-100 percentage
  chr.col = 1L, start.col = 2L, end.col = 3L,
  strand.col = 4L, coverage.col = 5L, freqC.col = 6L
)

# methylKit workflow: methRead -> unite -> calculateDiffMeth; returns NULL on failure.
run_methylkit <- function(sample_files, sample_ids, treatment, name) {
  tryCatch({
    mrl <- methylKit::methRead(
      location = as.list(sample_files),
      sample.id = as.list(sample_ids),
      assembly = opt$assembly,
      pipeline = pipeline_spec,
      header = TRUE, sep = "\t",
      context = opt[["methylkit-context"]],
      resolution = "base",
      treatment = treatment,
      mincov = 0
    )
    meth <- methylKit::unite(mrl, destrand = FALSE, min.per.group = NULL, mc.cores = opt[["mc-cores"]])
    diff <- methylKit::calculateDiffMeth(
      meth, overdispersion = "MN", adjust = opt[["methylkit-adjust"]],
      test = opt[["methylkit-test"]], mc.cores = opt[["mc-cores"]]
    )
    getData(diff)
  }, error = function(e) {
    message(sprintf("[%s] methylKit failed: %s", name, conditionMessage(e)))
    NULL
  })
}

# Single-genotype main pipeline: reads mock/treat, calculates effect size for each site,
# performs methylKit tests only on sites meeting per-replicate depth thresholds, returns one-sided p-values.
call_line_methylkit <- function(name, wt_rows, wt_ratio_mat, wt_depth_mat, mock_files, treat_files) {
  target_keys <- wt_rows$site_id
  allow_dup <- opt[["allow-duplicate-sites"]]

  mock_tables  <- lapply(mock_files,  read_adjusted, target_keys = target_keys, allow_dup = allow_dup)
  treat_tables <- lapply(treat_files, read_adjusted, target_keys = target_keys, allow_dup = allow_dup)

  # Complete-case: Site must be present in all mock and all treat files
  testable <- target_keys
  for (tb in mock_tables)  testable <- intersect(testable, tb$site_id)
  for (tb in treat_tables) testable <- intersect(testable, tb$site_id)

  if (length(testable) == 0) {
    return(list(name = name, results = list()))
  }

  w_idx <- match(testable, wt_rows$site_id)

  # Mutant Mock/Treat matrices. (Values for non-covered sites are 0.)
  mock_depth_mat  <- do.call(cbind, lapply(mock_tables,  function(tb) tb$depth[match(testable, tb$site_id)]))
  treat_ratio_mat <- do.call(cbind, lapply(treat_tables, function(tb) tb$ratio[match(testable, tb$site_id)]))
  treat_depth_mat <- do.call(cbind, lapply(treat_tables, function(tb) tb$depth[match(testable, tb$site_id)]))

  mut_ratio <- rowMeans(treat_ratio_mat)
  wt_ratio  <- wt_rows$wt_ratio[w_idx]
  delta     <- wt_ratio - mut_ratio

  # Depth threshold per replicate: WT/Mutant replicates must have depth >= min-depth;
  # any failing replicate excluded (status="depth_filter").
  min_d <- opt[["min-depth"]]
  depth_ok <- function(m) (rowSums(m >= min_d) == ncol(m)) & complete.cases(m)
  pass <- depth_ok(mock_depth_mat) & depth_ok(treat_depth_mat) &
          depth_ok(wt_depth_mat[w_idx, , drop = FALSE])

  # One-sided p-value, q-value (qvalue is double-sided FDR, adjusted by --methylkit-adjust);
  # sites not passing filters left as NA.
  p_one <- rep(NA_real_, length(testable))
  q_two <- rep(NA_real_, length(testable))
  mdiff_all <- rep(NA_real_, length(testable))
  names(p_one) <- names(q_two) <- names(mdiff_all) <- testable

  testable_pass <- testable[pass]
  if (length(testable_pass) > 0) {
    sel <- w_idx[pass]
    sites_g <- wt_rows[sel]
    wt_y <- wt_ratio_mat[sel, , drop = FALSE]
    wt_n <- wt_depth_mat[sel, , drop = FALSE]

    sample_dir <- file.path(opt$outdir, paste0(name, ".methylkit.samples"))
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    sample_files <- character(0)
    sample_ids   <- character(0)
    treatment    <- integer(0)
    for (i in seq_len(ncol(wt_y))) {
      fp <- file.path(sample_dir, sprintf("WT.rep%d.methyl.tsv", i))
      write_methyl_sample(fp, sites_g, wt_n[, i], wt_y[, i])
      sample_files <- c(sample_files, fp)
      sample_ids   <- c(sample_ids, sprintf("WT.rep%d", i))
      treatment    <- c(treatment, 1L)
    }
    for (j in seq_along(treat_tables)) {
      tb  <- treat_tables[[j]]
      idx <- match(sites_g$site_id, tb$site_id)
      fp  <- file.path(sample_dir, sprintf("%s.rep%d.methyl.tsv", name, j))
      write_methyl_sample(fp, sites_g, tb$depth[idx], tb$ratio[idx])
      sample_files <- c(sample_files, fp)
      sample_ids   <- c(sample_ids, sprintf("%s.rep%d", name, j))
      treatment    <- c(treatment, 0L)
    }

    test_res <- run_methylkit(sample_files, sample_ids, treatment, name)
    if (!is.null(test_res) && nrow(test_res) > 0) {
      res_key <- make_key(as.data.table(test_res), c("chr", "start", "strand"))
      hit <- match(testable_pass, res_key)
      ok <- !is.na(hit)
      p_two <- rep(NA_real_, length(testable_pass))
      mdiff <- rep(NA_real_, length(testable_pass))
      q_vec <- rep(NA_real_, length(testable_pass))
      p_two[ok]  <- test_res$pvalue[hit[ok]]
      mdiff[ok]  <- test_res$meth.diff[hit[ok]]
      q_vec[ok]  <- test_res$qvalue[hit[ok]]
      p_one[testable_pass] <- ifelse(!is.na(mdiff) & mdiff > 0, p_two / 2, 1)
      q_two[testable_pass] <- q_vec
      mdiff_all[testable_pass] <- mdiff
    }
  }

  use_q <- identical(opt[["significance"]], "q")
  results <- lapply(seq_along(testable), function(i) {
    k <- testable[i]
    status <- if (!pass[i]) {
      "depth_filter"
    } else if (is.na(p_one[[k]])) {
      "failed"
    } else {
      "ok"
    }
    # q-criterion: qvalue < alpha and positive direction; p-criterion uses the raw one-sided p.
    is_dep <- if (use_q) {
      pass[i] && !is.na(q_two[[k]]) && q_two[[k]] < opt$alpha &&
        !is.na(mdiff_all[[k]]) && mdiff_all[[k]] > 0 &&
        !is.na(delta[i]) && delta[i] >= opt[["min-delta"]]
    } else {
      pass[i] && !is.na(p_one[[k]]) && p_one[[k]] < opt$alpha &&
        !is.na(delta[i]) && delta[i] >= opt[["min-delta"]]
    }
    list(
      site_id = k, context = wt_rows$context[w_idx[i]], WT_Ratio = wt_ratio[i],
      Mock_Mean_Depth = mean(mock_depth_mat[i, ]), Mut_Ratio = mut_ratio[i], delta = delta[i],
      p = p_one[[k]], q = q_two[[k]], meth_diff = mdiff_all[[k]],
      status = status, is_dependent = is_dep
    )
  })

  list(name = name, results = results)
}

# ---------------------------------------------------------------------------
# Output-prefix inference (aligned with call_dependent.py / call_dependent_betabinomial.R)
# ---------------------------------------------------------------------------
infer_output_prefix <- function(g1, g2) {
  m1 <- sub("^(.+?)[_-]([12])([-_.].*)?$", "\\1\\3", g1)
  m2 <- sub("^(.+?)[_-]([12])([-_.].*)?$", "\\1\\3", g2)
  if (m1 == m2) return(m1)
  c1 <- strsplit(g1, "")[[1]]
  c2 <- strsplit(g2, "")[[1]]
  n <- min(length(c1), length(c2))
  i <- 1
  while (i <= n && c1[i] == c2[i]) i <- i + 1
  prefix <- paste0(c1[seq_len(max(0, i - 1))], collapse = "")
  prefix <- sub("[_\\-.]+$", "", prefix)
  if (nzchar(prefix)) return(prefix)
  paste0(g1, "_vs_", g2)
}

# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
main <- function() {
  suppressPackageStartupMessages(library(methylKit))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  dt_wt <- fread(opt$tsv, showProgress = FALSE)
  for (col in c(tsv_key_cols, opt[["ratio-col"]])) {
    if (!col %in% names(dt_wt)) stop(sprintf("--tsv missing column: %s", col))
  }

  wt_ratio_cols <- if (is.null(opt[["wt-ratio-cols"]]) || !nzchar(opt[["wt-ratio-cols"]])) {
    grep("^Treat-.*_ratio$", names(dt_wt), value = TRUE)
  } else {
    split_csv(opt[["wt-ratio-cols"]])
  }
  wt_depth_cols <- if (is.null(opt[["wt-depth-cols"]]) || !nzchar(opt[["wt-depth-cols"]])) {
    grep("^Treat-.*_depth$", names(dt_wt), value = TRUE)
  } else {
    split_csv(opt[["wt-depth-cols"]])
  }
  ratio_prefix <- sub("_ratio$", "", wt_ratio_cols)
  depth_prefix <- sub("_depth$", "", wt_depth_cols)
  common <- intersect(ratio_prefix, depth_prefix)
  if (length(common) == 0) {
    stop("No paired WT per-replicate Treat ratio/depth columns found in --tsv")
  }
  wt_ratio_cols <- paste0(common, "_ratio")
  wt_depth_cols <- paste0(common, "_depth")

  wt_rows <- data.table(
    site_id  = make_key(dt_wt, tsv_key_cols),
    context  = if (opt[["context-col"]] %in% names(dt_wt)) dt_wt[[opt[["context-col"]]]] else "NA",
    wt_ratio = as.numeric(dt_wt[[opt[["ratio-col"]]]])
  )
  # Keep the raw key-column values so write_methyl_sample can emit
  # chr/start/end/strand.
  for (j in seq_along(key_cols)) {
    wt_rows[[key_cols[j]]] <- dt_wt[[tsv_key_cols[j]]]
  }

  wt_ratio_mat <- matrix(NA_real_, nrow = nrow(dt_wt), ncol = length(common))
  wt_depth_mat <- matrix(NA_real_, nrow = nrow(dt_wt), ncol = length(common))
  for (i in seq_along(common)) {
    wt_ratio_mat[, i] <- as.numeric(dt_wt[[wt_ratio_cols[i]]])
    wt_depth_mat[, i] <- as.numeric(dt_wt[[wt_depth_cols[i]]])
  }

  line1 <- call_line_methylkit(
    opt$genotype1, wt_rows, wt_ratio_mat, wt_depth_mat,
    split_csv(opt[["genotype1-mock-files"]]),
    split_csv(opt[["genotype1-treat-files"]])
  )
  line2 <- call_line_methylkit(
    opt$genotype2, wt_rows, wt_ratio_mat, wt_depth_mat,
    split_csv(opt[["genotype2-mock-files"]]),
    split_csv(opt[["genotype2-treat-files"]])
  )

  output_prefix <- if (is.null(opt[["output-prefix"]]) || !nzchar(opt[["output-prefix"]])) {
    infer_output_prefix(opt$genotype1, opt$genotype2)
  } else {
    opt[["output-prefix"]]
  }

  p_label <- "p.value"
  q_label <- "q.value"

  # Line-level table (all sites for each mutant line)
  line_to_dt <- function(res) {
    rows <- rbindlist(lapply(res$results, function(r) data.table(
      site_id = r$site_id, context = r$context, WT_Ratio = r$WT_Ratio,
      Mock_Mean_Depth = r$Mock_Mean_Depth, Mut_Ratio = r$Mut_Ratio,
      delta = r$delta, p = r$p, q = r$q, meth_diff = r$meth_diff,
      status = r$status, is_dependent = r$is_dependent
    )))
    if (nrow(rows) == 0) return(rows)
    setnames(rows,
             c("delta", "p", "q", "meth_diff"),
             c("Delta_Ratio(WT-Mut)", p_label, q_label, "meth.diff"))
    keys_split <- tstrsplit(rows$site_id, "|", fixed = TRUE)
    for (j in seq_along(key_cols)) rows[[key_cols[j]]] <- keys_split[[j]]
    rows[, c(key_cols, setdiff(names(rows), c("site_id", key_cols))), with = FALSE]
  }

  fwrite(line_to_dt(line1), file.path(opt$outdir, paste0(opt$genotype1, ".dependent-sites.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(line_to_dt(line2), file.path(opt$outdir, paste0(opt$genotype2, ".dependent-sites.tsv")),
         sep = "\t", quote = FALSE)

  # Merge the two mutant lines (AND logic)
  l1_map <- setNames(lapply(line1$results, `[[`, "is_dependent"), sapply(line1$results, `[[`, "site_id"))
  l2_map <- setNames(lapply(line2$results, `[[`, "is_dependent"), sapply(line2$results, `[[`, "site_id"))
  l1_res <- setNames(line1$results, sapply(line1$results, `[[`, "site_id"))
  l2_res <- setNames(line2$results, sapply(line2$results, `[[`, "site_id"))
  common_keys <- sort(intersect(names(l1_map), names(l2_map)))

  final <- rbindlist(lapply(common_keys, function(k) {
    if (!(isTRUE(l1_map[[k]]) && isTRUE(l2_map[[k]]))) return(NULL)
    w <- wt_rows[site_id == k]
    mut_ratio <- (l1_res[[k]]$Mut_Ratio + l2_res[[k]]$Mut_Ratio) / 2
    data.table(
      site_id = k, context = w$context, WT_Ratio = w$wt_ratio,
      Mut_Ratio = mut_ratio, Delta = w$wt_ratio - mut_ratio,
      line1_p = l1_res[[k]]$p, line1_q = l1_res[[k]]$q,
      line2_p = l2_res[[k]]$p, line2_q = l2_res[[k]]$q
    )
  }))

  final_out <- file.path(opt$outdir, paste0(output_prefix, ".dependent.tsv"))
  final_cols <- c(key_cols, "context", "WT_Ratio", "Mut_Ratio", "Delta_Ratio(WT-Mut)",
                  paste0(opt$genotype1, "_", p_label), paste0(opt$genotype1, "_", q_label),
                  paste0(opt$genotype2, "_", p_label), paste0(opt$genotype2, "_", q_label))
  if (nrow(final) == 0) {
    empty <- as.data.table(setNames(
      rep(list(character(0)), length(final_cols)),
      final_cols
    ))
    fwrite(empty, final_out, sep = "\t", quote = FALSE)
    final_count <- 0
  } else {
    setnames(final, "Delta", "Delta_Ratio(WT-Mut)")
    setnames(final, "line1_p", paste0(opt$genotype1, "_", p_label))
    setnames(final, "line1_q", paste0(opt$genotype1, "_", q_label))
    setnames(final, "line2_p", paste0(opt$genotype2, "_", p_label))
    setnames(final, "line2_q", paste0(opt$genotype2, "_", q_label))
    keys_split <- tstrsplit(final$site_id, "|", fixed = TRUE)
    for (j in seq_along(key_cols)) final[[key_cols[j]]] <- keys_split[[j]]
    final <- final[, c(key_cols, "context", "WT_Ratio", "Mut_Ratio",
                       "Delta_Ratio(WT-Mut)",
                       paste0(opt$genotype1, "_", p_label), paste0(opt$genotype1, "_", q_label),
                       paste0(opt$genotype2, "_", p_label), paste0(opt$genotype2, "_", q_label)),
                   with = FALSE]
    fwrite(final, final_out, sep = "\t", quote = FALSE)
    final_count <- nrow(final)
  }

  summary <- list(
    tsv = opt$tsv,
    key_cols = key_cols,
    tsv_key_cols = tsv_key_cols,
    ratio_col = opt[["ratio-col"]],
    adjusted_ratio_col = opt[["adjusted-ratio-col"]],
    adjusted_depth_col = opt[["adjusted-depth-col"]],
    min_depth = opt[["min-depth"]],
    min_delta = opt[["min-delta"]],
    method = paste0("methylKit_", opt[["methylkit-test"]]),
    methylkit_adjust = opt[["methylkit-adjust"]],
    significance = opt[["significance"]],
    alpha = opt$alpha,
    methylkit_context = opt[["methylkit-context"]],
    wt_candidate_count = nrow(wt_rows),
    wt_replicate_count = length(common),
    genotype1 = opt$genotype1,
    genotype1_all_sites_count = length(line1$results),
    genotype1_dependent_count = sum(sapply(line1$results, `[[`, "is_dependent")),
    genotype2 = opt$genotype2,
    genotype2_all_sites_count = length(line2$results),
    genotype2_dependent_count = sum(sapply(line2$results, `[[`, "is_dependent")),
    final_all_sites_count = length(common_keys),
    final_dependent_count = final_count
  )
  summary_out <- file.path(opt$outdir, paste0(output_prefix, ".dependent.summary.json"))
  writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), summary_out)

  message(sprintf("WT candidates: %d", nrow(wt_rows)))
  message(sprintf("%s all sites: %d", opt$genotype1, length(line1$results)))
  message(sprintf("%s dependent sites: %d", opt$genotype1, sum(sapply(line1$results, `[[`, "is_dependent"))))
  message(sprintf("%s all sites: %d", opt$genotype2, length(line2$results)))
  message(sprintf("%s dependent sites: %d", opt$genotype2, sum(sapply(line2$results, `[[`, "is_dependent"))))
  message(sprintf("Final dependent sites: %d", final_count))
  message(sprintf("Final output: %s", final_out))
}

main()
