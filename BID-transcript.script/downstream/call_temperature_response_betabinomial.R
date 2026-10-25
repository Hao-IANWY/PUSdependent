#!/usr/bin/env Rscript
# ============================================================================
# call_temperature_response_betabinomial.R — beta-binomial LRT 版温度响应位点判定
#
# 与 call_temperature_response_methylKit.R（methylKit 版）并列：比较同一 WT
# 材料在两个温度下的修饰差异（group2=WT_27 vs group1=WT_21），找出修饰比例
# 随温度变化的位点。输入接口、候选位点口径（两温度 fisher 阳性并集，由
# run.sh 预合并传入 --tsv）、每重复深度门槛、效应量口径与判定三条件与其
# 完全一致，仅把统计检验换成 beta-binomial（方法与 call_dependent_betabinomial.R
# 同源，温度响应不预设方向，故用双侧 P，方向单独记录在 direction 列）：
#
#   1) 全局离散度 theta：在全部过深度门的位点上拟合含 site 随机截距的
#      beta-binomial GLMM（cbind(gap, depth-gap) ~ group + (1|site)），
#      收敛/Hessian 异常时用 BFGS 重试，仍失败则全部待测位点记
#      status="theta_failed"；
#   2) 逐位点检验：固定全局 theta（map 锁定 betadisp）在同一 data 上拟合
#      cbind(gap, depth-gap) ~ group 与 ~ 1 两个模型，--site-test 决定双侧
#      P 的来源：lrt（默认）取似然比 D = 2(ll_full − ll_null) 的 χ²(1) 检验，
#      wald 取 full 模型 group 系数的 Wald z（分离位点 SE 爆炸、P 退化，
#      仅复现用）。β/SE/OR 审计列两种口径下均取自 full 模型 Wald；
#   3) BH 校正家族 = 全部候选位点：深度不达标的行 P 为 NA，不参与排序
#      但计入 n（与 call_dependent_betabinomial.R 一致，偏保守）；
#   4) 边界处理：两组全部重复均为 0%（或均为 100%）且 delta=0 的无对比
#      位点保守取 P=1（status="boundary_p1"）；"一组非 0、另一组全 0"的
#      完全分离位点无特判逻辑：lrt 口径下检验有效，wald 口径下 P 会退化。
#
# 判定（三个条件同时满足才判 temperature-responsive）：
#   1. 每重复深度门槛：两个温度组的 Mock 与 Treat 每一个重复 depth 都
#      >= --min-depth（默认 10）；
#   2. 显著性：默认 q = 双侧 P 的 BH FDR < --alpha；--significance p 可退回
#      原始双侧 p 口径；
#   3. |Delta_Ratio(group2-group1)| >= --min-delta（fraction）。
#   group ratio = 该温度组 Treat 重复 ratio 的简单均值；direction 列单独
#   记录方向（up_in_<group2> / down_in_<group2>），不参与判定。
#
# 输出口径与 call_temperature_response_methylKit.R 完全同构：all-sites.tsv /
# response.tsv / summary.json 同名同列（q.value 仍在 direction 之后，下游
# data_prep 按固定列号读 direction），行级表额外追加 beta_coef / beta_se /
# odds_ratio 三列审计信息，theta 写入 summary JSON。
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
  library(parallel)
})

# glmmTMB 延迟加载：脚本主体用到时再 require。
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
              help = "Per-site two-sided P source: lrt (default; likelihood-ratio test of group, valid under complete separation) or wald (legacy Wald z; degenerates on separation, kept for reproducing methylKit-era behavior comparisons)."),
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
# 读取 adjusted TSV，只保留目标位点的 (ratio, depth, gap)
# ---------------------------------------------------------------------------
read_adjusted <- function(path, target_keys, allow_dup) {
  # adjusted 文件可能高达数千万行；先用 grep -F 按 (chr,pos,strand) 精确预过滤，
  # 只把目标位点所在行交给 fread，避免全量读取。（key 中的 "\t" 天然锚定了 pos/strand
  # 边界，chr 为唯一转录本 ID，不存在跨位点子串误匹配。）
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
# 统计核心（与 call_dependent_betabinomial.R 同源，genotype 换成 group）
# ---------------------------------------------------------------------------
validate_counts <- function(gap, depth, label) {
  bad <- !is.na(gap) & !is.na(depth) &
    (gap < 0 | depth < 0 | gap > depth | gap != round(gap) | depth != round(depth))
  if (any(bad)) {
    stop(sprintf("%s: %d invalid gap/depth count(s)", label, sum(bad)))
  }
}

# 全局离散度：包含 site 随机截距的 beta-binomial GLMM；失败返回 NA
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
  # 收敛或 Hessian 异常时改用 BFGS 重试
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

# 逐位点拟合：固定全局 theta，返回双侧 P（--site-test 决定 Wald/LRT 口径）、
# 系数与 odds ratio
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
  # 每组至少需要 2 个有效观测
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

  # 同一 data、同一锁死的 theta 下拟合给定公式；失败用 BFGS 重试
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

  # 双侧 P：lrt 取 2(ll_full − ll_null) 的 χ²(1) 检验，完全分离位点统计量
  # 有限、检验有效；wald 取 full 模型 Wald z，分离位点 SE 爆炸、P 退化
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
# 主流程
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

  # complete-case：位点需同时出现在所有 mock 与所有 treat 文件（对齐 dependent 版）
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

  # 每重复深度门槛：两个温度组的 Mock 与 Treat 每一个重复 depth 都须
  # >= min-depth；任何一个重复不达标即剔除（status="depth_filter"）。
  min_d <- opt[["min-depth"]]
  depth_ok <- function(m) (rowSums(m >= min_d) == ncol(m)) & complete.cases(m)
  pass <- depth_ok(fetch_cols(g1_mock, "depth")) & depth_ok(fetch_cols(g1_treat, "depth")) &
          depth_ok(fetch_cols(g2_mock, "depth")) & depth_ok(fetch_cols(g2_treat, "depth"))

  # Treat ratio/depth/gap 矩阵。adjusted 全量网格对未覆盖位点写 depth=0 且
  # ratio/gap=0（防除零约定），已由上面的每重复门槛整体挡在检验之外。
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

  # 双侧 P 值、BH q 值（家族 = 全部候选位点，深度不达标行为 NA 但计入 n）、
  # 审计用系数；深度不达标或检验未运行的位点留 NA
  n_testable <- length(testable)
  p_two <- rep(NA_real_, n_testable)
  q_two <- rep(NA_real_, n_testable)
  beta_v <- rep(NA_real_, n_testable)
  se_v   <- rep(NA_real_, n_testable)
  or_v   <- rep(NA_real_, n_testable)
  status_extra <- rep(NA_character_, n_testable)  # pass 位点的细分状态
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

      # 边界处理：两组全部重复均为 0%（或均为 100%）且 delta=0 的无对比位点
      # 保守取 P=1、β=0、OR=1（只覆盖拟合失败的位点）
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

  # BH 家族 = 全部候选位点（NA 不参与排序但计入 n，见头注）
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
    # q 口径：双侧 P 的 BH FDR < alpha；p 口径：原始双侧 P < alpha。
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

  # q.value 放在 direction 之后：下游 data_prep 仍按固定列号读 direction
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
