#!/usr/bin/env Rscript
# ============================================================================
# call_dependent_betabinomial.R — beta-binomial GLMM 版 dependent 位点判定
#
# 与 call_dependent_methylKit.R（methylKit 版）并列：输入接口、深度门槛、
# 效应量口径、两株系 AND 判定与输出表结构与其完全一致，仅把统计检验换成
# Jie Huang（Beacon）dependency_stats_standalone.R 的方法（见仓库
# dependent_jiehuang/extracted/dependency_stats_standalone/scripts/）：
#
#   1) 全局离散度 theta：在全部过深度门的位点上拟合含 site 随机截距的
#      beta-binomial GLMM（cbind(gap, depth-gap) ~ genotype + (1|site)），
#      收敛/Hessian 异常时用 BFGS 重试，仍失败则该株系全部待测位点记
#      status="theta_failed"（原方法此时整体报错退出，这里保持逐株系降级）；
#   2) 逐位点检验：固定全局 theta（map 锁定 betadisp）在同一 data 上拟合
#      cbind(gap, depth-gap) ~ genotype 与 ~ 1 两个模型，--site-test 决定双侧
#      P 的来源：lrt（默认）取似然比 D = 2(ll_full − ll_null) 的 χ²(1) 检验，
#      wald 取 full 模型 genotypemutant 系数的 Wald z（旧口径）。两者都按
#      full 模型 β 的符号折单侧 loss P：β<0（突变体修饰下降方向）取双侧
#      P/2，否则取 1−P/2。LRT 在完全分离（如突变体三重复全 0%）下统计量
#      有限、检验有效；Wald 在分离位点 SE 爆炸、P 退化为 ~0.5（仅复现旧
#      结果用）。β/SE/OR 审计列两种口径下均取自 full 模型 Wald；
#   3) BH 校正家族 = 该株系全部候选位点：深度不达标的行 P 为 NA，不参与
#      排序但计入 n，因此比只在校验池内校正更保守（与原方法一致）；
#   4) 边界处理：两组全部重复均为 0%（或均为 100%）且 delta=0 的无对比
#      位点保守取 P=1（status="boundary_p1"）；"WT 非 0、突变体全 0"的
#      完全分离位点无特判逻辑：lrt 口径下似然比统计量有限、检验有效，
#      wald 口径下 P 会退化（原方法的已知敏感点）。
#
# 与原方法的两处有意差异（保持本仓库 dependent 判定框架不变，保证与
# methylKit 版逐位点可比）：
#   - 深度门槛沿用本流程严格口径：突变株 Mock 与 Treat 的每一个重复、WT
#     每个重复都须 depth >= --min-depth（原方法只要求每基因型 >=2 个文库
#     达标）；Mock 仍只用于深度门槛与 Mock_Mean_Depth 列，不进入模型；
#   - delta 沿用均值差口径（WT_Ratio − Mut_Ratio，见 methylKit 版头注），
#     原方法用 pooled 差；两者只在个别边界情形不同。
# 本流程的输入表没有批次注释，模型不含 batch 项（原方法支持可选 batch）。
#
# 显著性口径（--significance）：
#   q（默认）  单侧 P 的 BH FDR < --alpha 且 delta >= --min-delta；
#   p          原始单侧 P < --alpha（仅用于复现/对照）。
# 输出口径与 call_dependent_methylKit.R 完全同构：p.value 列为单侧 loss P，
# q.value 为其 BH 校正；行级表额外给出 β、SE、odds ratio 供审计，theta
# 写入 summary JSON。
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
  make_option(c("--wt-gap-cols"), type = "character", default = NULL,
              help = "Comma-separated WT per-replicate Treat gap columns. Default: auto-detect '^Treat-.*_gap$'."),
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
  make_option(c("--adjusted-gap-col"), type = "character", default = "gap",
              help = "Gap (modified-read count) column in adjusted TSV files. Default: gap."),
  make_option(c("--min-depth"), type = "double", default = 10,
              help = "Minimum per-replicate depth required for every Mock/Treat (and WT) replicate. Default: 10."),
  make_option(c("--min-delta"), type = "double", default = 0.05,
              help = "Effect-size floor (WT_Ratio - Mut_Ratio), in fraction. Default: 0.05."),
  make_option(c("--alpha"), type = "double", default = 0.05,
              help = "Significance cutoff for --significance. Default: 0.05."),
  make_option(c("--significance"), type = "character", default = "q",
              help = "Statistic used for dependent calls: q = BH FDR of one-sided P (default), p = raw one-sided P."),
  make_option(c("--site-test"), type = "character", default = "lrt",
              help = "Per-site two-sided P source: lrt (default; likelihood-ratio test of genotype, valid under complete separation) or wald (legacy Wald z; degenerates on separation, kept for reproducing old outputs)."),
  make_option(c("--cores"), type = "integer", default = 1,
              help = "Cores for the per-site beta-binomial fits (mclapply). Default: 1."),
  make_option(c("--allow-duplicate-sites"), action = "store_true", default = FALSE,
              help = "Keep first occurrence of duplicate keys instead of failing.")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("tsv", "genotype1", "genotype1-treat-files", "genotype1-mock-files",
              "genotype2", "genotype2-treat-files", "genotype2-mock-files", "outdir")
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
tsv_key_cols <- if (is.null(opt[["tsv-key-cols"]]) || !nzchar(opt[["tsv-key-cols"]])) key_cols else split_csv(opt[["tsv-key-cols"]])
if (length(tsv_key_cols) != length(key_cols)) {
  stop("--tsv-key-cols and --key-cols must define the same number of columns")
}

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
# 统计核心（移植自 dependency_stats_standalone.R，去掉本流程用不到的 batch 项）
# ---------------------------------------------------------------------------
validate_counts <- function(gap, depth, label) {
  bad <- !is.na(gap) & !is.na(depth) &
    (gap < 0 | depth < 0 | gap > depth | gap != round(gap) | depth != round(depth))
  if (any(bad)) {
    stop(sprintf("%s: %d invalid gap/depth count(s)", label, sum(bad)))
  }
}

# 单侧 loss P：β<0 取双侧 P/2，否则取 1 - 双侧 P/2
loss_one_sided_p <- function(beta, two_sided_p) {
  ifelse(beta < 0, two_sided_p / 2, 1 - two_sided_p / 2)
}

# 全局离散度：包含 site 随机截距的 beta-binomial GLMM；失败返回 NA
estimate_theta <- function(site_ids, wg, wd, mg, md) {
  n_site <- length(site_ids)
  n_wt_lib <- ncol(wg)
  long <- data.frame(
    site = factor(rep(site_ids, n_wt_lib + ncol(mg))),
    gap = c(as.vector(wg), as.vector(mg)),
    depth = c(as.vector(wd), as.vector(md)),
    genotype = factor(
      c(rep("WT", n_wt_lib * n_site), rep("mutant", ncol(mg) * n_site)),
      levels = c("WT", "mutant")
    )
  )
  long <- long[!is.na(long$gap) & !is.na(long$depth) & long$depth > 0, , drop = FALSE]

  formula <- cbind(gap, depth - gap) ~ genotype + (1 | site)
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

# 逐位点拟合：固定全局 theta，返回单侧 P（--site-test 决定 Wald/LRT 口径）、
# 系数与 odds ratio
fit_site <- function(wg, wd, mg, md, theta) {
  data <- data.frame(
    gap = c(wg, mg), depth = c(wd, md),
    genotype = factor(
      c(rep("WT", length(wg)), rep("mutant", length(mg))),
      levels = c("WT", "mutant")
    )
  )
  data <- droplevels(data[!is.na(data$gap) & !is.na(data$depth) &
                            data$depth > 0, , drop = FALSE])

  empty <- c(
    raw_p = NA_real_, two_sided_p = NA_real_, coefficient = NA_real_,
    standard_error = NA_real_, odds_ratio = NA_real_, fit_ok = 0
  )
  # 每个基因型至少需要 2 个有效观测
  if (sum(data$genotype == "WT") < 2L || sum(data$genotype == "mutant") < 2L) {
    return(empty)
  }
  full_formula <- cbind(gap, depth - gap) ~ genotype
  null_formula <- cbind(gap, depth - gap) ~ 1
  design <- model.matrix(full_formula, data)
  if (qr(design)$rank < ncol(design)) {
    return(empty)
  }

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
  if (!"genotypemutant" %in% rownames(coefficients)) {
    return(empty)
  }
  estimate <- coefficients["genotypemutant", "Estimate"]
  se <- coefficients["genotypemutant", "Std. Error"]
  if (!all(is.finite(c(estimate, se)))) {
    return(empty)
  }

  # 双侧 P：lrt 取 2(ll_full − ll_null) 的 χ²(1) 检验，完全分离位点统计量
  # 有限、检验有效；wald 取 full 模型 Wald z，分离位点 SE 爆炸、P 退化
  two_sided_p <- if (identical(opt[["site-test"]], "wald")) {
    coefficients["genotypemutant", "Pr(>|z|)"]
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
    raw_p = loss_one_sided_p(estimate, two_sided_p),
    two_sided_p = two_sided_p,
    coefficient = estimate,
    standard_error = se,
    odds_ratio = exp(estimate),
    fit_ok = 1
  )
}

# ---------------------------------------------------------------------------
# 单株系主流程：读 mock/treat，逐位点算效应量，对每重复深度达标的位点估计
# 全局 theta 后逐位点拟合 beta-binomial，返回单侧 P 与 BH q 值。
# ---------------------------------------------------------------------------
call_line_betabinomial <- function(name, wt_rows, wt_ratio_mat, wt_depth_mat, wt_gap_mat,
                                   mock_files, treat_files) {
  target_keys <- wt_rows$site_id
  allow_dup <- opt[["allow-duplicate-sites"]]

  mock_tables  <- lapply(mock_files,  read_adjusted, target_keys = target_keys, allow_dup = allow_dup)
  treat_tables <- lapply(treat_files, read_adjusted, target_keys = target_keys, allow_dup = allow_dup)

  # complete-case：位点需同时出现在所有 mock 与所有 treat 文件（对齐 methylKit 版）
  testable <- target_keys
  for (tb in mock_tables)  testable <- intersect(testable, tb$site_id)
  for (tb in treat_tables) testable <- intersect(testable, tb$site_id)

  if (length(testable) == 0) {
    return(list(name = name, results = list(), theta = NA_real_))
  }

  w_idx <- match(testable, wt_rows$site_id)

  # 突变株 Mock/Treat 的 ratio/depth/gap 矩阵。adjusted 是全量位点网格，未覆盖位点
  # depth=0 且 ratio/gap 硬编码为 0，由下面的每重复深度门槛整体挡在检验之外。
  mock_depth_mat  <- do.call(cbind, lapply(mock_tables,  function(tb) tb$depth[match(testable, tb$site_id)]))
  treat_ratio_mat <- do.call(cbind, lapply(treat_tables, function(tb) tb$ratio[match(testable, tb$site_id)]))
  treat_depth_mat <- do.call(cbind, lapply(treat_tables, function(tb) tb$depth[match(testable, tb$site_id)]))
  treat_gap_mat   <- do.call(cbind, lapply(treat_tables, function(tb) tb$gap[match(testable, tb$site_id)]))
  validate_counts(treat_gap_mat, treat_depth_mat, sprintf("%s Treat", name))

  mut_ratio <- rowMeans(treat_ratio_mat)
  wt_ratio  <- wt_rows$wt_ratio[w_idx]
  delta     <- wt_ratio - mut_ratio

  # 每重复深度门槛：突变株 Mock 与 Treat 的全部重复、WT 的全部重复都须
  # >= min-depth；任何一个重复不达标即剔除（status="depth_filter"）。
  min_d <- opt[["min-depth"]]
  depth_ok <- function(m) (rowSums(m >= min_d) == ncol(m)) & complete.cases(m)
  pass <- depth_ok(mock_depth_mat) & depth_ok(treat_depth_mat) &
          depth_ok(wt_depth_mat[w_idx, , drop = FALSE])

  # 单侧 P 值、BH q 值（家族 = 全部候选位点，深度不达标行为 NA 但计入 n）、
  # 审计用系数；深度不达标或检验未运行的位点留 NA
  n_testable <- length(testable)
  p_one <- rep(NA_real_, n_testable)
  q_val <- rep(NA_real_, n_testable)
  beta_v <- rep(NA_real_, n_testable)
  se_v   <- rep(NA_real_, n_testable)
  or_v   <- rep(NA_real_, n_testable)
  status_extra <- rep(NA_character_, n_testable)  # pass 位点的细分状态
  names(p_one) <- names(q_val) <- names(beta_v) <- names(se_v) <-
    names(or_v) <- names(status_extra) <- testable

  theta <- NA_real_
  testable_pass <- testable[pass]
  if (length(testable_pass) > 0) {
    sel <- which(pass)
    wg <- wt_gap_mat[w_idx[sel], , drop = FALSE]
    wd <- wt_depth_mat[w_idx[sel], , drop = FALSE]
    mg <- treat_gap_mat[sel, , drop = FALSE]
    md <- treat_depth_mat[sel, , drop = FALSE]
    wr <- wt_ratio_mat[w_idx[sel], , drop = FALSE]
    mr <- treat_ratio_mat[sel, , drop = FALSE]
    delta_pass <- delta[sel]

    theta <- estimate_theta(testable_pass, wg, wd, mg, md)
    if (!is.finite(theta) || theta <= 0) {
      message(sprintf("[%s] global theta estimation failed; all tested sites marked theta_failed", name))
      status_extra[sel] <- "theta_failed"
    } else {
      fits <- parallel::mclapply(seq_len(nrow(wg)), function(i) {
        fit_site(wg[i, ], wd[i, ], mg[i, ], md[i, ], theta)
      }, mc.cores = opt[["cores"]], mc.preschedule = TRUE)
      fit_mat <- do.call(rbind, fits)

      # 边界处理：两组全部重复均为 0%（或均为 100%）且 delta=0 的无对比位点
      # 保守取 P=1、β=0、OR=1（与原方法一致；只覆盖拟合失败的位点）
      observation_n <- ncol(wr) + ncol(mr)
      all_one <- rowSums(!is.na(wr) & wr == 1) + rowSums(!is.na(mr) & mr == 1) ==
        observation_n
      all_zero <- rowSums(!is.na(wr) & wr == 0) + rowSums(!is.na(mr) & mr == 0) ==
        observation_n
      boundary <- fit_mat[, "fit_ok"] != 1 & (all_one | all_zero) &
        !is.na(delta_pass) & delta_pass == 0
      fit_mat[boundary, "raw_p"] <- 1
      fit_mat[boundary, "two_sided_p"] <- 1
      fit_mat[boundary, "coefficient"] <- 0
      fit_mat[boundary, "odds_ratio"] <- 1

      p_one[sel] <- fit_mat[, "raw_p"]
      beta_v[sel] <- fit_mat[, "coefficient"]
      se_v[sel] <- fit_mat[, "standard_error"]
      or_v[sel] <- fit_mat[, "odds_ratio"]
      status_extra[sel] <- ifelse(boundary, "boundary_p1",
                                  ifelse(fit_mat[, "fit_ok"] == 1, "ok", "fit_failed"))
    }
  }

  # BH 家族 = 全部候选位点（NA 不参与排序但计入 n，见头注）
  q_val <- stats::p.adjust(p_one, method = "BH", n = n_testable)

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
    # q 口径：单侧 P 的 BH FDR < alpha；p 口径：原始单侧 P < alpha。
    # 方向已内嵌在单侧 P 中（β<0 才可能得到小 P；boundary P=1 自然落选）。
    is_dep <- if (use_q) {
      pass[i] && !is.na(q_val[[k]]) && q_val[[k]] < opt$alpha &&
        !is.na(delta[i]) && delta[i] >= opt[["min-delta"]]
    } else {
      pass[i] && !is.na(p_one[[k]]) && p_one[[k]] < opt$alpha &&
        !is.na(delta[i]) && delta[i] >= opt[["min-delta"]]
    }
    list(
      site_id = k, context = wt_rows$context[w_idx[i]], WT_Ratio = wt_ratio[i],
      Mock_Mean_Depth = mean(mock_depth_mat[i, ]), Mut_Ratio = mut_ratio[i], delta = delta[i],
      p = p_one[[k]], q = q_val[[k]], beta = beta_v[[k]], se = se_v[[k]],
      or = or_v[[k]], status = status, is_dependent = is_dep
    )
  })

  list(name = name, results = results, theta = theta)
}

# ---------------------------------------------------------------------------
# 输出前缀推断（对齐 call_dependent.py / call_dependent_methylKit.R）
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
# 主流程
# ---------------------------------------------------------------------------
main <- function() {
  suppressPackageStartupMessages(library(glmmTMB))
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
  wt_gap_cols <- if (is.null(opt[["wt-gap-cols"]]) || !nzchar(opt[["wt-gap-cols"]])) {
    grep("^Treat-.*_gap$", names(dt_wt), value = TRUE)
  } else {
    split_csv(opt[["wt-gap-cols"]])
  }
  ratio_prefix <- sub("_ratio$", "", wt_ratio_cols)
  depth_prefix <- sub("_depth$", "", wt_depth_cols)
  gap_prefix <- sub("_gap$", "", wt_gap_cols)
  common <- Reduce(intersect, list(ratio_prefix, depth_prefix, gap_prefix))
  if (length(common) == 0) {
    stop("No paired WT per-replicate Treat ratio/depth/gap columns found in --tsv")
  }
  wt_ratio_cols <- paste0(common, "_ratio")
  wt_depth_cols <- paste0(common, "_depth")
  wt_gap_cols <- paste0(common, "_gap")

  wt_rows <- data.table(
    site_id  = make_key(dt_wt, tsv_key_cols),
    context  = if (opt[["context-col"]] %in% names(dt_wt)) dt_wt[[opt[["context-col"]]]] else "NA",
    wt_ratio = as.numeric(dt_wt[[opt[["ratio-col"]]]])
  )

  n_wt_rep <- length(common)
  wt_ratio_mat <- matrix(NA_real_, nrow = nrow(dt_wt), ncol = n_wt_rep)
  wt_depth_mat <- matrix(NA_real_, nrow = nrow(dt_wt), ncol = n_wt_rep)
  wt_gap_mat <- matrix(NA_real_, nrow = nrow(dt_wt), ncol = n_wt_rep)
  for (i in seq_len(n_wt_rep)) {
    wt_ratio_mat[, i] <- as.numeric(dt_wt[[wt_ratio_cols[i]]])
    wt_depth_mat[, i] <- as.numeric(dt_wt[[wt_depth_cols[i]]])
    wt_gap_mat[, i] <- as.numeric(dt_wt[[wt_gap_cols[i]]])
  }
  validate_counts(wt_gap_mat, wt_depth_mat, "WT Treat")

  line1 <- call_line_betabinomial(
    opt$genotype1, wt_rows, wt_ratio_mat, wt_depth_mat, wt_gap_mat,
    split_csv(opt[["genotype1-mock-files"]]),
    split_csv(opt[["genotype1-treat-files"]])
  )
  line2 <- call_line_betabinomial(
    opt$genotype2, wt_rows, wt_ratio_mat, wt_depth_mat, wt_gap_mat,
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

  # line-level 表（每个株系所有位点）
  line_to_dt <- function(res) {
    rows <- rbindlist(lapply(res$results, function(r) data.table(
      site_id = r$site_id, context = r$context, WT_Ratio = r$WT_Ratio,
      Mock_Mean_Depth = r$Mock_Mean_Depth, Mut_Ratio = r$Mut_Ratio,
      delta = r$delta, p = r$p, q = r$q, beta = r$beta, se = r$se,
      or = r$or, status = r$status, is_dependent = r$is_dependent
    )))
    if (nrow(rows) == 0) return(rows)
    setnames(rows,
             c("delta", "p", "q", "beta", "se", "or"),
             c("Delta_Ratio(WT-Mut)", p_label, q_label, "beta_coef", "beta_se", "odds_ratio"))
    keys_split <- tstrsplit(rows$site_id, "|", fixed = TRUE)
    for (j in seq_along(key_cols)) rows[[key_cols[j]]] <- keys_split[[j]]
    rows[, c(key_cols, setdiff(names(rows), c("site_id", key_cols))), with = FALSE]
  }

  fwrite(line_to_dt(line1), file.path(opt$outdir, paste0(opt$genotype1, ".dependent-sites.tsv")),
         sep = "\t", quote = FALSE)
  fwrite(line_to_dt(line2), file.path(opt$outdir, paste0(opt$genotype2, ".dependent-sites.tsv")),
         sep = "\t", quote = FALSE)

  # 合并两个株系（AND）
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

  count_status <- function(res, level) sum(sapply(res$results, function(r) r$status == level))
  summary <- list(
    tsv = opt$tsv,
    key_cols = key_cols,
    tsv_key_cols = tsv_key_cols,
    ratio_col = opt[["ratio-col"]],
    adjusted_ratio_col = opt[["adjusted-ratio-col"]],
    adjusted_depth_col = opt[["adjusted-depth-col"]],
    adjusted_gap_col = opt[["adjusted-gap-col"]],
    min_depth = opt[["min-depth"]],
    min_delta = opt[["min-delta"]],
    method = "betaBinomial_glmmTMB",
    significance = opt[["significance"]],
    site_test = opt[["site-test"]],
    alpha = opt$alpha,
    wt_candidate_count = nrow(wt_rows),
    wt_replicate_count = n_wt_rep,
    genotype1 = opt$genotype1,
    genotype1_theta = unname(line1$theta),
    genotype1_all_sites_count = length(line1$results),
    genotype1_dependent_count = sum(sapply(line1$results, `[[`, "is_dependent")),
    genotype1_fit_failed_count = count_status(line1, "fit_failed"),
    genotype1_boundary_p1_count = count_status(line1, "boundary_p1"),
    genotype1_theta_failed = count_status(line1, "theta_failed") > 0,
    genotype2 = opt$genotype2,
    genotype2_theta = unname(line2$theta),
    genotype2_all_sites_count = length(line2$results),
    genotype2_dependent_count = sum(sapply(line2$results, `[[`, "is_dependent")),
    genotype2_fit_failed_count = count_status(line2, "fit_failed"),
    genotype2_boundary_p1_count = count_status(line2, "boundary_p1"),
    genotype2_theta_failed = count_status(line2, "theta_failed") > 0,
    final_all_sites_count = length(common_keys),
    final_dependent_count = final_count
  )
  summary_out <- file.path(opt$outdir, paste0(output_prefix, ".dependent.summary.json"))
  writeLines(toJSON(summary, pretty = TRUE, auto_unbox = TRUE), summary_out)

  message(sprintf("WT candidates: %d", nrow(wt_rows)))
  for (ln in list(line1, line2)) {
    message(sprintf("%s theta: %s", ln$name,
                    if (is.finite(ln$theta)) sprintf("%.4f", ln$theta) else "NA"))
    message(sprintf("%s all sites: %d", ln$name, length(ln$results)))
    message(sprintf("%s dependent sites: %d", ln$name,
                    sum(sapply(ln$results, `[[`, "is_dependent"))))
  }
  message(sprintf("Final dependent sites: %d", final_count))
  message(sprintf("Final output: %s", final_out))
}

main()
