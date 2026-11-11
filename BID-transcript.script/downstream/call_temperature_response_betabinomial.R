#!/usr/bin/env Rscript
# ============================================================================
# call_temperature_response_betabinomial.R - temperature-responsive site identification via beta-binomial LRT
#
# Companion to call_temperature_response_methylKit.R. Identifies sites where modification ratio 
# changes significantly between two temperatures (group2=WT_27 vs group1=WT_21).
# Statistical framework uses beta-binomial (GLMM) following methods in call_dependent_betabinomial.R.
# Temperature response testing uses two-sided P-values, with direction stored separately.
#
# 1) Global dispersion (theta): Fit beta-binomial GLMM (cbind(gap, depth-gap) ~ group + (1|site))
#    on sites meeting depth filters. Fallback to BFGS on convergence failure.
# 2) Per-site testing: Using locked global theta, fit models (group vs intercept).
#    --site-test selects P-value source: 'lrt' (default; chi-squared(1) likelihood-ratio test, robust)
#    or 'wald' (Wald z-statistic, for comparison runs).
# 3) Multiple testing correction: BH-FDR applied across all candidate sites.
# 4) Boundary conditions: Sites with 0% or 100% modification in both groups and delta=0 
#    conservatively assigned P=1.
#
# Identification criteria (all must be met):
#   1. Depth threshold: Every Mock/Treat replicate in both temperature groups >= --min-depth (default 10).
#   2. Significance: BH FDR of two-sided P < --alpha (default).
#   3. Effect size: |Delta_Ratio(group2-group1)| >= --min-delta.
#
# Output:
#   all-sites.tsv, response.tsv, summary.json (same structure as methylKit version).
#   Includes beta_coef, beta_se, odds_ratio for audit; summary includes theta.
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(parallel)
})

# Load glmmTMB lazily; the main workflow attaches it when needed.
if (!requireNamespace("glmmTMB", quietly = TRUE)) {
  stop("This script requires the 'glmmTMB' package.")
}

option_list <- list(
  make_option(c("--tsv"), type = "character",
              help = "Candidate site TSV with a header and coordinate + context columns (union of both temperatures' positives)."),
  make_option(c("--key-cols"), type = "character", default = "chr,pos,strand",
              help = "Comma-separated site key columns in --tsv and adjusted TSV files. Default: chr,pos,strand."),
  make_option(c("--context-col"), type = "character", default = "context",
              help = "Optional context column name. Default: context."),
  make_option(c("--group1"), type = "character",
              help = "Name of the baseline group, e.g. WT_21."),
  make_option(c("--group1-treat-files"), type = "character",
              help = "Comma-separated Treat adjusted TSV files for group1."),
  make_option(c("--group1-mock-files"), type = "character",
              help = "Comma-separated Mock adjusted TSV files for group1."),
  make_option(c("--group2"), type = "character",
              help = "Name of the comparison group, e.g. WT_27."),
  make_option(c("--group2-treat-files"), type = "character",
              help = "Comma-separated Treat adjusted TSV files for group2."),
  make_option(c("--group2-mock-files"), type = "character",
              help = "Comma-separated Mock adjusted TSV files for group2."),
  make_option(c("--outdir"), type = "character",
              help = "Output directory for all-sites, response, and summary files."),
  make_option(c("--output-prefix"), type = "character", default = NULL,
              help = "Prefix for output files. Default: '<group2>_vs_<group1>'."),
  make_option(c("--adjusted-ratio-col"), type = "character", default = "ratio",
              help = "Ratio column in adjusted TSV files. Default: ratio."),
  make_option(c("--adjusted-depth-col"), type = "character", default = "depth",
              help = "Depth column in adjusted TSV files. Default: depth."),
  make_option(c("--adjusted-gap-col"), type = "character", default = "gap",
              help = "Gap (modified-read count) column in adjusted TSV files. Default: gap."),
  make_option(c("--min-depth"), type = "double", default = 10,
              help = "Minimum per-replicate depth required for every Mock/Treat replicate of both groups. Default: 10."),
  make_option(c("--min-delta"), type = "double", default = 0.05,
              help = "Effect-size floor |ratio(group2) - ratio(group1)|, in fraction. Default: 0.05."),
  make_option(c("--alpha"), type = "double", default = 0.05,
              help = "Significance cutoff for --significance. Default: 0.05."),
  make_option(c("--significance"), type = "character", default = "q",
              help = "Statistic used for response calls: q = BH FDR of two-sided P (default), p = raw two-sided P."),
  make_option(c("--site-test"), type = "character", default = "lrt",
              help = "Per-site two-sided P source: lrt (default; likelihood-ratio test of group, valid under complete separation) or wald (Wald z; degenerates on separation, retained for reproducing wald-based results)."),
  make_option(c("--cores"), type = "integer", default = 1,
              help = "Cores for the per-site beta-binomial fits (mclapply). Default: 1."),
  make_option(c("--allow-duplicate-sites"), action = "store_true", default = FALSE,
              help = "Keep first occurrence of duplicate keys instead of failing.")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("tsv", "group1", "group1-treat-files", "group1-mock-files",
              "group2", "group2-treat-files", "group2-mock-files", "outdir")
for (r in required) {
  if (is.null(opt[[r]]) || !nzchar(opt[[r]])) stop(sprintf("Missing required option: --%s", r))
}
if (!opt[["significance"]] %in% c("p", "q")) {
  stop("--significance must be 'p' or 'q'")
}
if (!opt[["site-test"]] %in% c("wald", "lrt")) {
  stop("--site-test must be 'wald' or 'lrt'")
}
if (!(opt[["min-delta"]] >= 0 && opt[["min-delta"]] <= 1)) {
  stop("--min-delta must be a fraction in [0, 1]")
}
if (is.na(opt[["cores"]]) || opt[["cores"]] < 1L) {
  stop("--cores must be a positive integer")
}

split_csv <- function(x) {
  items <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  items[nzchar(items)]
}

key_cols <- split_csv(opt[["key-cols"]])

make_key <- function(dt, cols) {
  do.call(paste, c(as.list(dt)[cols], list(sep = "|")))
}

# ---------------------------------------------------------------------------
# Read an adjusted TSV, keeping only (ratio, depth, gap) for target sites
# ---------------------------------------------------------------------------
read_adjusted <- function(path, target_keys, allow_dup) {
  # Adjusted files can reach tens of millions of rows, so pre-filter exactly
  # with grep -F on (chr,pos,strand) and let fread see only target-site rows
  # instead of the whole file. ("\t" inside each key anchors the pos/strand
  # boundaries, and chr is a unique transcript ID, so cross-site substring
  # matches cannot occur.)
  key_file <- tempfile(fileext = ".keys")
  writeLines(gsub("|", "\t", target_keys, fixed = TRUE), key_file)
  hdr <- strsplit(readLines(path, n = 1), "\t")[[1]]
  cmd <- sprintf("LC_ALL=C grep -F -f %s %s", shQuote(key_file), shQuote(path))
  dt <- fread(cmd = cmd, header = FALSE, col.names = hdr, showProgress = FALSE)
  setnames(dt, opt[["adjusted-ratio-col"]], "ratio")
  setnames(dt, opt[["adjusted-depth-col"]], "depth")
  setnames(dt, opt[["adjusted-gap-col"]], "gap")
  dt[, site_id := make_key(.SD, key_cols)]
  dt <- dt[site_id %in% target_keys]
  dt[, ratio := as.numeric(ratio)]
  dt[, depth := as.numeric(depth)]
  dt[, gap := as.numeric(gap)]
  if (!allow_dup && anyDuplicated(dt$site_id)) {
    stop(sprintf("%s: duplicate site key", path))
  }
  dt[, .(site_id, ratio, depth, gap)]
}

# ---------------------------------------------------------------------------
# Statistical core (shared with call_dependent_betabinomial.R, with genotype
# replaced by group)
# ---------------------------------------------------------------------------
validate_counts <- function(gap, depth, label) {
  bad <- !is.na(gap) & !is.na(depth) &
    (gap < 0 | depth < 0 | gap > depth | gap != round(gap) | depth != round(depth))
  if (any(bad)) {
    stop(sprintf("%s: %d invalid gap/depth count(s)", label, sum(bad)))
  }
}

# Global dispersion: beta-binomial GLMM with a site random intercept;
# returns NA on failure.
estimate_theta <- function(site_ids, g1g, g1d, g2g, g2d) {
  n_site <- length(site_ids)
  n1 <- ncol(g1g)
  long <- data.frame(
    site = factor(rep(site_ids, n1 + ncol(g2g))),
    gap = c(as.vector(g1g), as.vector(g2g)),
    depth = c(as.vector(g1d), as.vector(g2d)),
    grp = factor(
      c(rep(opt$group1, n1 * n_site), rep(opt$group2, ncol(g2g) * n_site)),
      levels = c(opt$group1, opt$group2)
    )
  )
  long <- long[!is.na(long$gap) & !is.na(long$depth) & long$depth > 0, , drop = FALSE]

  formula <- cbind(gap, depth - gap) ~ grp + (1 | site)
  fit <- tryCatch(
    suppressWarnings(glmmTMB::glmmTMB(
      formula, family = glmmTMB::betabinomial(), data = long
    )),
    error = function(e) NULL
  )
  # Retry with BFGS on convergence or Hessian failure.
  if (is.null(fit) || fit$fit$convergence != 0L || !isTRUE(fit$sdr$pdHess)) {
    fit <- tryCatch(
      suppressWarnings(glmmTMB::glmmTMB(
        formula, family = glmmTMB::betabinomial(), data = long,
        control = glmmTMB::glmmTMBControl(
          optimizer = stats::optim, optArgs = list(method = "BFGS")
        )
      )),
      error = function(e) NULL
    )
  }
  if (is.null(fit) || fit$fit$convergence != 0L || !isTRUE(fit$sdr$pdHess)) {
    return(NA_real_)
  }
  as.numeric(stats::sigma(fit))
}

# Per-site fit with the global theta locked; returns the two-sided P
# (--site-test selects the Wald or LRT criteria), coefficient, and odds ratio.
fit_site <- function(g1g, g1d, g2g, g2d, theta) {
  data <- data.frame(
    gap = c(g1g, g2g), depth = c(g1d, g2d),
    grp = factor(
      c(rep(opt$group1, length(g1g)), rep(opt$group2, length(g2g))),
      levels = c(opt$group1, opt$group2)
    )
  )
  data <- droplevels(data[!is.na(data$gap) & !is.na(data$depth) &
                            data$depth > 0, , drop = FALSE])

  empty <- c(
    two_sided_p = NA_real_, coefficient = NA_real_,
    standard_error = NA_real_, odds_ratio = NA_real_, fit_ok = 0
  )
  # Each group needs at least 2 valid observations.
  if (sum(data$grp == opt$group1) < 2L || sum(data$grp == opt$group2) < 2L) {
    return(empty)
  }
  full_formula <- cbind(gap, depth - gap) ~ grp
  null_formula <- cbind(gap, depth - gap) ~ 1
  design <- model.matrix(full_formula, data)
  if (qr(design)$rank < ncol(design)) {
    return(empty)
  }
  coef_name <- paste0("grp", opt$group2)

  # Fit the given formula on the same data with theta locked; retry with
  # BFGS on failure.
  do_fit <- function(formula, control = glmmTMB::glmmTMBControl()) {
    tryCatch(
      suppressWarnings(glmmTMB::glmmTMB(
        formula, family = glmmTMB::betabinomial(), data = data,
        start = list(betadisp = log(theta)),
        map = list(betadisp = factor(NA)),
        control = control
      )),
      error = function(e) NULL
    )
  }
  fit_or_retry <- function(formula) {
    fit <- do_fit(formula)
    if (is.null(fit) || fit$fit$convergence != 0L || !isTRUE(fit$sdr$pdHess)) {
      fit <- do_fit(formula, glmmTMB::glmmTMBControl(
        optimizer = stats::optim, optArgs = list(method = "BFGS")
      ))
    }
    fit
  }
  fit_full <- fit_or_retry(full_formula)
  if (is.null(fit_full) || fit_full$fit$convergence != 0L ||
        !isTRUE(fit_full$sdr$pdHess)) {
    return(empty)
  }

  coefficients <- summary(fit_full)$coefficients$cond
  if (!coef_name %in% rownames(coefficients)) {
    return(empty)
  }
  estimate <- coefficients[coef_name, "Estimate"]
  se <- coefficients[coef_name, "Std. Error"]
  if (!all(is.finite(c(estimate, se)))) {
    return(empty)
  }

  # Two-sided P: 'lrt' tests 2(ll_full - ll_null) against chi-squared(1); the
  # statistic stays finite under complete separation, so the test remains
  # valid. 'wald' uses the full-model Wald z, whose SE explodes at separated
  # sites and degenerates the P-value.
  two_sided_p <- if (identical(opt[["site-test"]], "wald")) {
    coefficients[coef_name, "Pr(>|z|)"]
  } else {
    fit_null <- fit_or_retry(null_formula)
    if (is.null(fit_null) || fit_null$fit$convergence != 0L ||
          !isTRUE(fit_null$sdr$pdHess)) {
      return(empty)
    }
    lrt_stat <- max(0, 2 * (as.numeric(stats::logLik(fit_full)) -
                              as.numeric(stats::logLik(fit_null))))
    stats::pchisq(lrt_stat, df = 1, lower.tail = FALSE)
  }
  if (!is.finite(two_sided_p)) {
    return(empty)
  }

  c(
    two_sided_p = two_sided_p,
    coefficient = estimate,
    standard_error = se,
    odds_ratio = exp(estimate),
    fit_ok = 1
  )
}

# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
main <- function() {
  suppressPackageStartupMessages(library(glmmTMB))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  dt_cand <- fread(opt$tsv, showProgress = FALSE)
  for (col in key_cols) {
    if (!col %in% names(dt_cand)) stop(sprintf("--tsv missing column: %s", col))
  }
  cand <- data.table(site_id = make_key(dt_cand, key_cols))
  cand[, context := if (opt[["context-col"]] %in% names(dt_cand)) dt_cand[[opt[["context-col"]]]] else "NA"]
  for (j in seq_along(key_cols)) cand[[key_cols[j]]] <- dt_cand[[key_cols[j]]]
  if (!opt[["allow-duplicate-sites"]] && anyDuplicated(cand$site_id)) {
    stop(sprintf("%s: duplicate site key", opt$tsv))
  }
  target_keys <- cand$site_id
  allow_dup <- opt[["allow-duplicate-sites"]]

  g1_treat <- lapply(split_csv(opt[["group1-treat-files"]]), read_adjusted, target_keys = target_keys, allow_dup = allow_dup)
  g1_mock  <- lapply(split_csv(opt[["group1-mock-files"]]),  read_adjusted, target_keys = target_keys, allow_dup = allow_dup)
  g2_treat <- lapply(split_csv(opt[["group2-treat-files"]]), read_adjusted, target_keys = target_keys, allow_dup = allow_dup)
  g2_mock  <- lapply(split_csv(opt[["group2-mock-files"]]),  read_adjusted, target_keys = target_keys, allow_dup = allow_dup)

  # Complete-case analysis: a site must appear in every mock and every treat
  # file (matching the dependent-site version).
  testable <- target_keys
  for (tb in c(g1_mock, g1_treat, g2_mock, g2_treat)) testable <- intersect(testable, tb$site_id)

  group1 <- opt$group1
  group2 <- opt$group2
  output_prefix <- if (is.null(opt[["output-prefix"]]) || !nzchar(opt[["output-prefix"]])) {
    paste0(group2, "_vs_", group1)
  } else {
    opt[["output-prefix"]]
  }
  all_out  <- file.path(opt$outdir, paste0(output_prefix, ".all-sites.tsv"))
  resp_out <- file.path(opt$outdir, paste0(output_prefix, ".response.tsv"))
  json_out <- file.path(opt$outdir, paste0(output_prefix, ".response.summary.json"))

  ratio1_col <- paste0(group1, "_Ratio")
  ratio2_col <- paste0(group2, "_Ratio")
  delta_col  <- sprintf("Delta_Ratio(%s-%s)", group2, group1)
  depth1_col <- paste0("Mock_Mean_Depth_", group1)
  depth2_col <- paste0("Mock_Mean_Depth_", group2)

  write_empty_outputs <- function() {
    empty_all <- as.data.table(setNames(
      rep(list(character(0)), length(c(key_cols, "context", ratio1_col, ratio2_col, delta_col,
                                       depth1_col, depth2_col, "p.value", "status", "is_responsive",
                                       "direction", "q.value", "beta_coef", "beta_se", "odds_ratio"))),
      c(key_cols, "context", ratio1_col, ratio2_col, delta_col,
        depth1_col, depth2_col, "p.value", "status", "is_responsive",
        "direction", "q.value", "beta_coef", "beta_se", "odds_ratio")
    ))
    fwrite(empty_all, all_out, sep = "\t", quote = FALSE)
    empty_resp <- empty_all[, c(key_cols, "context", ratio1_col, ratio2_col, delta_col,
                                "p.value", "direction", "q.value"), with = FALSE]
    fwrite(empty_resp, resp_out, sep = "\t", quote = FALSE)
  }

  if (length(testable) == 0) {
    write_empty_outputs()
    summary <- list(tsv = opt$tsv, key_cols = key_cols, group1 = group1, group2 = group2,
                    candidate_count = nrow(cand), testable_count = 0L, responsive_count = 0L,
                    up_count = 0L, down_count = 0L, theta = NA)
    writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), json_out)
    message("No testable sites (empty complete-case intersection).")
    return(invisible(NULL))
  }

  c_idx <- match(testable, cand$site_id)

  fetch_cols <- function(tables, field) {
    do.call(cbind, lapply(tables, function(tb) tb[[field]][match(testable, tb$site_id)]))
  }

  # Per-replicate depth filter: every Mock and Treat replicate of both
  # temperature groups must have depth >= min-depth; a site failing any
  # replicate is dropped (status="depth_filter").
  min_d <- opt[["min-depth"]]
  depth_ok <- function(m) (rowSums(m >= min_d) == ncol(m)) & complete.cases(m)
  pass <- depth_ok(fetch_cols(g1_mock, "depth")) & depth_ok(fetch_cols(g1_treat, "depth")) &
          depth_ok(fetch_cols(g2_mock, "depth")) & depth_ok(fetch_cols(g2_treat, "depth"))

  # Treat ratio/depth/gap matrices. The adjusted full-site grid writes
  # depth = 0 and ratio/gap = 0 for uncovered sites (a division-guard
  # convention); the per-replicate filter above already excludes them.
  tr1_y <- fetch_cols(g1_treat, "ratio")
  tr1_n <- fetch_cols(g1_treat, "depth")
  tr1_g <- fetch_cols(g1_treat, "gap")
  tr2_y <- fetch_cols(g2_treat, "ratio")
  tr2_n <- fetch_cols(g2_treat, "depth")
  tr2_g <- fetch_cols(g2_treat, "gap")
  validate_counts(tr1_g, tr1_n, sprintf("%s Treat", group1))
  validate_counts(tr2_g, tr2_n, sprintf("%s Treat", group2))

  ratio1 <- rowMeans(tr1_y)
  ratio2 <- rowMeans(tr2_y)
  delta  <- ratio2 - ratio1

  mock_depth1 <- rowMeans(fetch_cols(g1_mock, "depth"))
  mock_depth2 <- rowMeans(fetch_cols(g2_mock, "depth"))

  # Two-sided P and BH q values (family = all candidate sites; depth-filtered
  # rows stay NA but count towards n) plus audit coefficients; sites failing
  # the depth filter or not tested remain NA.
  n_testable <- length(testable)
  p_two <- rep(NA_real_, n_testable)
  q_two <- rep(NA_real_, n_testable)
  beta_v <- rep(NA_real_, n_testable)
  se_v   <- rep(NA_real_, n_testable)
  or_v   <- rep(NA_real_, n_testable)
  status_extra <- rep(NA_character_, n_testable)  # fine-grained status for passing sites
  names(p_two) <- names(q_two) <- names(beta_v) <- names(se_v) <-
    names(or_v) <- names(status_extra) <- testable

  theta <- NA_real_
  testable_pass <- testable[pass]
  if (length(testable_pass) > 0) {
    sel <- which(pass)
    g1g <- tr1_g[sel, , drop = FALSE]
    g1d <- tr1_n[sel, , drop = FALSE]
    g2g <- tr2_g[sel, , drop = FALSE]
    g2d <- tr2_n[sel, , drop = FALSE]
    r1 <- tr1_y[sel, , drop = FALSE]
    r2 <- tr2_y[sel, , drop = FALSE]
    delta_pass <- delta[sel]

    theta <- estimate_theta(testable_pass, g1g, g1d, g2g, g2d)
    if (!is.finite(theta) || theta <= 0) {
      message(sprintf("[%s] global theta estimation failed; all tested sites marked theta_failed", output_prefix))
      status_extra[sel] <- "theta_failed"
    } else {
      fits <- parallel::mclapply(seq_len(nrow(g1g)), function(i) {
        fit_site(g1g[i, ], g1d[i, ], g2g[i, ], g2d[i, ], theta)
      }, mc.cores = opt[["cores"]], mc.preschedule = TRUE)
      fit_mat <- do.call(rbind, fits)

      # Boundary handling: sites without contrast (all replicates 0%, or all
      # 100%, in both groups) and delta = 0 are conservatively assigned P = 1,
      # beta = 0, and OR = 1 (only fit-failed sites are covered).
      observation_n <- ncol(r1) + ncol(r2)
      all_one <- rowSums(!is.na(r1) & r1 == 1) + rowSums(!is.na(r2) & r2 == 1) ==
        observation_n
      all_zero <- rowSums(!is.na(r1) & r1 == 0) + rowSums(!is.na(r2) & r2 == 0) ==
        observation_n
      boundary <- fit_mat[, "fit_ok"] != 1 & (all_one | all_zero) &
        !is.na(delta_pass) & delta_pass == 0
      fit_mat[boundary, "two_sided_p"] <- 1
      fit_mat[boundary, "coefficient"] <- 0
      fit_mat[boundary, "odds_ratio"] <- 1

      p_two[sel] <- fit_mat[, "two_sided_p"]
      beta_v[sel] <- fit_mat[, "coefficient"]
      se_v[sel] <- fit_mat[, "standard_error"]
      or_v[sel] <- fit_mat[, "odds_ratio"]
      status_extra[sel] <- ifelse(boundary, "boundary_p1",
                                  ifelse(fit_mat[, "fit_ok"] == 1, "ok", "fit_failed"))
    }
  }

  # BH family = all candidate sites
  q_two <- stats::p.adjust(p_two, method = "BH", n = n_testable)

  use_q <- identical(opt[["significance"]], "q")
  results <- lapply(seq_len(n_testable), function(i) {
    k <- testable[i]
    status <- if (!pass[i]) {
      "depth_filter"
    } else if (!is.na(status_extra[[k]])) {
      status_extra[[k]]
    } else {
      "failed"
    }
    direction <- if (is.na(delta[i])) {
      NA_character_
    } else if (delta[i] > 0) {
      paste0("up_in_", group2)
    } else if (delta[i] < 0) {
      paste0("down_in_", group2)
    } else {
      "unchanged"
    }
    # q-criterion: two-sided P BH FDR < alpha; p-criterion: raw two-sided P < alpha.
    is_resp <- if (use_q) {
      pass[i] && !is.na(q_two[[k]]) && q_two[[k]] < opt$alpha &&
        !is.na(delta[i]) && abs(delta[i]) >= opt[["min-delta"]]
    } else {
      pass[i] && !is.na(p_two[[k]]) && p_two[[k]] < opt$alpha &&
        !is.na(delta[i]) && abs(delta[i]) >= opt[["min-delta"]]
    }
    list(
      site_id = k, context = cand$context[c_idx[i]],
      ratio1 = ratio1[i], ratio2 = ratio2[i], delta = delta[i],
      depth1 = mock_depth1[i], depth2 = mock_depth2[i],
      p = p_two[[k]], q = q_two[[k]], beta = beta_v[[k]], se = se_v[[k]],
      or = or_v[[k]], status = status, is_responsive = is_resp,
      direction = if (is_resp) direction else NA_character_
    )
  })

  all_dt <- rbindlist(lapply(results, function(r) data.table(
    site_id = r$site_id, context = r$context, ratio1 = r$ratio1, ratio2 = r$ratio2,
    delta = r$delta, depth1 = r$depth1, depth2 = r$depth2,
    p = r$p, status = r$status, is_responsive = r$is_responsive,
    direction = r$direction, q = r$q, beta = r$beta, se = r$se, or = r$or
  )))
  setnames(all_dt, c("ratio1", "ratio2", "delta", "depth1", "depth2", "p", "q",
                     "beta", "se", "or"),
           c(ratio1_col, ratio2_col, delta_col, depth1_col, depth2_col,
             "p.value", "q.value", "beta_coef", "beta_se", "odds_ratio"))
  keys_split <- tstrsplit(all_dt$site_id, "|", fixed = TRUE)
  for (j in seq_along(key_cols)) all_dt[[key_cols[j]]] <- keys_split[[j]]
  all_dt <- all_dt[, c(key_cols, setdiff(names(all_dt), c("site_id", key_cols))), with = FALSE]
  fwrite(all_dt, all_out, sep = "\t", quote = FALSE)

  # q.value follows direction because downstream data_prep still reads the
  # direction column by fixed index.
  resp_dt <- all_dt[is_responsive == TRUE, c(key_cols, "context", ratio1_col, ratio2_col,
                                             delta_col, "p.value", "direction", "q.value"),
                    with = FALSE]
  fwrite(resp_dt, resp_out, sep = "\t", quote = FALSE)

  up_count   <- sum(resp_dt$direction == paste0("up_in_", group2))
  down_count <- sum(resp_dt$direction == paste0("down_in_", group2))

  count_status <- function(level) sum(sapply(results, function(r) r$status == level))
  summary <- list(
    tsv = opt$tsv,
    key_cols = key_cols,
    adjusted_ratio_col = opt[["adjusted-ratio-col"]],
    adjusted_depth_col = opt[["adjusted-depth-col"]],
    adjusted_gap_col = opt[["adjusted-gap-col"]],
    min_depth = opt[["min-depth"]],
    min_delta = opt[["min-delta"]],
    method = "betaBinomial_glmmTMB_two-sided",
    significance = opt[["significance"]],
    site_test = opt[["site-test"]],
    alpha = opt$alpha,
    candidate_count = nrow(cand),
    testable_count = n_testable,
    depth_filtered_count = sum(!pass),
    tested_count = count_status("ok"),
    theta_failed_count = count_status("theta_failed"),
    fit_failed_count = count_status("fit_failed"),
    boundary_p1_count = count_status("boundary_p1"),
    theta = unname(theta),
    group1 = group1,
    group1_replicate_count = length(g1_treat),
    group2 = group2,
    group2_replicate_count = length(g2_treat),
    responsive_count = nrow(resp_dt),
    up_count = up_count,
    down_count = down_count
  )
  writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), json_out)

  message(sprintf("Candidates: %d", nrow(cand)))
  message(sprintf("Testable sites: %d", n_testable))
  message(sprintf("Depth filtered: %d", sum(!pass)))
  message(sprintf("Theta: %s", if (is.finite(theta)) sprintf("%.4f", theta) else "NA"))
  message(sprintf("Responsive sites: %d (up_in_%s: %d, down_in_%s: %d)",
                  nrow(resp_dt), group2, up_count, group2, down_count))
  message(sprintf("Response output: %s", resp_out))
}

main()
