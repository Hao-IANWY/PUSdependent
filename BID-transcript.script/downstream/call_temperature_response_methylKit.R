#!/usr/bin/env Rscript
# ============================================================================
# call_temperature_response_methylKit.R — methylKit 版温度响应位点判定
#
# 与 call_dependent_methylKit.R（WT vs 突变株，单侧“修饰丢失”）并列，改为
# 比较同一 WT 材料在两个温度下的修饰差异（group2=WT_27 vs group1=WT_21），
# 找出修饰比例随温度发生变化的位点。检验流程同样对齐 Akalin et al. 2012：
#   methRead() -> unite() -> calculateDiffMeth()
#   多重复走 logReg（coverage 为权的二项 logistic 回归 + MN 过离散 F 检验）；
#   仅当总共 2 个样本（每组 1 个）时 methylKit 自动强制 fast.fisher。
#
# 分组约定：group2 记 treatment=1、group1 记 treatment=0，因此
#     meth.diff = (group2 加权平均修饰) - (group1 加权平均修饰)。
# 与 dependent 版不同，这里不预设方向——温度响应可以是升高或降低，故直接
# 使用 methylKit 的双侧 p 值，方向单独记录在 direction 列。
#
# 候选位点：两个温度 fisher 阳性位点的并集（仅单一温度阳性的位点正是
# “另一温度修饰丢失”的强候选），由 run.sh 预先合并传入 --tsv。
#
# 判定（三个条件同时满足才判 temperature-responsive）：
#   1. 每重复深度门槛：两个温度组的 Mock 与 Treat 每一个重复 depth 都
#      >= --min-depth（默认 10；排除“某温度下转录本不表达”造成的假修饰
#      差异，同时保证零/低覆盖的假观测不进检验）；
#   2. 双侧 p < --alpha；
#   3. |Delta_Ratio(group2-group1)| >= --min-delta（fraction）。
#   group ratio = 该温度组 Treat 重复 ratio 的简单均值。
#
# adjusted 是全量位点网格，未覆盖位点写 depth=0 且 ratio=0（adjustGap 防除
# 零约定，非有效观测）；上述每重复门槛保证这类伪观测被整体挡在检验之外，
# 因此无需零覆盖特判逻辑。
# ============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(jsonlite)
})

# methylKit 延迟加载：脚本主体用到时再 require，避免甲基化库缺失时连 help 都看不到。
if (!requireNamespace("methylKit", quietly = TRUE)) {
  stop("This script requires the 'methylKit' package. Install with:\n",
       "  conda install -c bioconda bioconductor-methylkit")
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
  make_option(c("--min-depth"), type = "double", default = 10,
              help = "Minimum per-replicate depth required for every Mock/Treat replicate of both groups. Default: 10."),
  make_option(c("--min-delta"), type = "double", default = 0.05,
              help = "Effect-size floor |ratio(group2) - ratio(group1)|, in fraction. Default: 0.05."),
  make_option(c("--alpha"), type = "double", default = 0.05,
              help = "Significance cutoff for the two-sided p-value. Default: 0.05."),
  make_option(c("--methylkit-test"), type = "character", default = "F",
              help = "methylKit test statistic for the multi-replicate logistic-regression path: F (default; overdispersion-corrected F test, degrades to Chisq when phi<=1) or Chisq. Only used when total samples > 2; with exactly 2 samples methylKit forces fast.fisher regardless."),
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

required <- c("tsv", "group1", "group1-treat-files", "group1-mock-files",
              "group2", "group2-treat-files", "group2-mock-files", "outdir")
for (r in required) {
  if (is.null(opt[[r]]) || !nzchar(opt[[r]])) stop(sprintf("Missing required option: --%s", r))
}
if (!opt[["methylkit-test"]] %in% c("F", "Chisq")) {
  stop("--methylkit-test must be 'F' or 'Chisq'")
}
if (!(opt[["min-delta"]] >= 0 && opt[["min-delta"]] <= 1)) {
  stop("--min-delta must be a fraction in [0, 1]")
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
# 读取 adjusted TSV，只保留目标位点的 (ratio, depth)
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
  dt[, site_id := make_key(.SD, key_cols)]
  dt <- dt[site_id %in% target_keys]
  dt[, ratio := as.numeric(ratio)]
  dt[, depth := as.numeric(depth)]
  if (!allow_dup && anyDuplicated(dt$site_id)) {
    stop(sprintf("%s: duplicate site key", path))
  }
  dt[, .(site_id, ratio, depth)]
}

# ---------------------------------------------------------------------------
# 把样本写成 methylKit 泛型格式（chr start end strand coverage freqC）
# freqC 用 fraction（0-1），methylKit 内部 numCs = round(coverage * freqC)
# ---------------------------------------------------------------------------
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

# methylKit 泛型读取的列映射
pipeline_spec <- list(
  fraction = TRUE,   # freqC 是 fraction（0-1），非 0-100 百分数
  chr.col = 1L, start.col = 2L, end.col = 3L,
  strand.col = 4L, coverage.col = 5L, freqC.col = 6L
)

# ---------------------------------------------------------------------------
# methylKit 工作流：methRead -> unite -> calculateDiffMeth；失败时返回 NULL。
# ---------------------------------------------------------------------------
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
      meth, overdispersion = "MN", adjust = "SLIM",
      test = opt[["methylkit-test"]], mc.cores = opt[["mc-cores"]]
    )
    getData(diff)
  }, error = function(e) {
    message(sprintf("[%s] methylKit failed: %s", name, conditionMessage(e)))
    NULL
  })
}

# ---------------------------------------------------------------------------
# 主流程
# ---------------------------------------------------------------------------
main <- function() {
  suppressPackageStartupMessages(library(methylKit))
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
                                       depth1_col, depth2_col, "p.value", "status", "is_responsive", "direction"))),
      c(key_cols, "context", ratio1_col, ratio2_col, delta_col,
        depth1_col, depth2_col, "p.value", "status", "is_responsive", "direction")
    ))
    fwrite(empty_all, all_out, sep = "\t", quote = FALSE)
    empty_resp <- empty_all[, c(key_cols, "context", ratio1_col, ratio2_col, delta_col, "p.value", "direction"), with = FALSE]
    fwrite(empty_resp, resp_out, sep = "\t", quote = FALSE)
  }

  if (length(testable) == 0) {
    write_empty_outputs()
    summary <- list(tsv = opt$tsv, key_cols = key_cols, group1 = group1, group2 = group2,
                    candidate_count = nrow(cand), testable_count = 0L, responsive_count = 0L)
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

  # Treat ratio/depth 矩阵。adjusted 全量网格对未覆盖位点写 depth=0 且
  # ratio=0（防除零约定），已由上面的每重复门槛整体挡在检验之外。
  tr1_n <- fetch_cols(g1_treat, "depth")
  tr1_y <- fetch_cols(g1_treat, "ratio")
  tr2_n <- fetch_cols(g2_treat, "depth")
  tr2_y <- fetch_cols(g2_treat, "ratio")

  ratio1 <- rowMeans(tr1_y)
  ratio2 <- rowMeans(tr2_y)
  delta  <- ratio2 - ratio1

  mock_depth1 <- rowMeans(fetch_cols(g1_mock, "depth"))
  mock_depth2 <- rowMeans(fetch_cols(g2_mock, "depth"))

  # 双侧 p 值（温度响应不预设方向）；深度不达标的位点留 NA
  p_two <- rep(NA_real_, length(testable))
  names(p_two) <- testable

  testable_pass <- testable[pass]
  if (length(testable_pass) > 0) {
    rows_pass <- which(pass)
    sites_g <- cand[c_idx[rows_pass]]
    g1_y_sub <- tr1_y[rows_pass, , drop = FALSE]
    g1_n_sub <- tr1_n[rows_pass, , drop = FALSE]
    g2_y_sub <- tr2_y[rows_pass, , drop = FALSE]
    g2_n_sub <- tr2_n[rows_pass, , drop = FALSE]

    sample_dir <- file.path(opt$outdir, paste0(output_prefix, ".methylkit.samples"))
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    sample_files <- character(0)
    sample_ids   <- character(0)
    treatment    <- integer(0)
    for (i in seq_len(ncol(g1_y_sub))) {
      fp <- file.path(sample_dir, sprintf("%s.rep%d.methyl.tsv", group1, i))
      write_methyl_sample(fp, sites_g, g1_n_sub[, i], g1_y_sub[, i])
      sample_files <- c(sample_files, fp)
      sample_ids   <- c(sample_ids, sprintf("%s.rep%d", group1, i))
      treatment    <- c(treatment, 0L)
    }
    for (j in seq_len(ncol(g2_y_sub))) {
      fp <- file.path(sample_dir, sprintf("%s.rep%d.methyl.tsv", group2, j))
      write_methyl_sample(fp, sites_g, g2_n_sub[, j], g2_y_sub[, j])
      sample_files <- c(sample_files, fp)
      sample_ids   <- c(sample_ids, sprintf("%s.rep%d", group2, j))
      treatment    <- c(treatment, 1L)
    }

    test_res <- run_methylkit(sample_files, sample_ids, treatment, output_prefix)
    if (!is.null(test_res) && nrow(test_res) > 0) {
      res_key <- make_key(as.data.table(test_res), c("chr", "start", "strand"))
      hit <- match(testable_pass, res_key)
      ok <- !is.na(hit)
      p_two[testable_pass[ok]] <- test_res$pvalue[hit[ok]]
    }
  }

  results <- lapply(seq_along(testable), function(i) {
    k <- testable[i]
    status <- if (!pass[i]) {
      "depth_filter"
    } else if (is.na(p_two[[k]])) {
      "failed"
    } else {
      "ok"
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
    is_resp <- pass[i] && !is.na(p_two[[k]]) && p_two[[k]] < opt$alpha &&
               !is.na(delta[i]) && abs(delta[i]) >= opt[["min-delta"]]
    list(
      site_id = k, context = cand$context[c_idx[i]],
      ratio1 = ratio1[i], ratio2 = ratio2[i], delta = delta[i],
      depth1 = mock_depth1[i], depth2 = mock_depth2[i],
      p = p_two[[k]], status = status, is_responsive = is_resp,
      direction = if (is_resp) direction else NA_character_
    )
  })

  all_dt <- rbindlist(lapply(results, function(r) data.table(
    site_id = r$site_id, context = r$context, ratio1 = r$ratio1, ratio2 = r$ratio2,
    delta = r$delta, depth1 = r$depth1, depth2 = r$depth2,
    p = r$p, status = r$status, is_responsive = r$is_responsive, direction = r$direction
  )))
  setnames(all_dt, c("ratio1", "ratio2", "delta", "depth1", "depth2", "p"),
           c(ratio1_col, ratio2_col, delta_col, depth1_col, depth2_col, "p.value"))
  keys_split <- tstrsplit(all_dt$site_id, "|", fixed = TRUE)
  for (j in seq_along(key_cols)) all_dt[[key_cols[j]]] <- keys_split[[j]]
  all_dt <- all_dt[, c(key_cols, setdiff(names(all_dt), c("site_id", key_cols))), with = FALSE]
  fwrite(all_dt, all_out, sep = "\t", quote = FALSE)

  resp_dt <- all_dt[is_responsive == TRUE, c(key_cols, "context", ratio1_col, ratio2_col,
                                             delta_col, "p.value", "direction"), with = FALSE]
  fwrite(resp_dt, resp_out, sep = "\t", quote = FALSE)

  up_count   <- sum(resp_dt$direction == paste0("up_in_", group2))
  down_count <- sum(resp_dt$direction == paste0("down_in_", group2))

  summary <- list(
    tsv = opt$tsv,
    key_cols = key_cols,
    adjusted_ratio_col = opt[["adjusted-ratio-col"]],
    adjusted_depth_col = opt[["adjusted-depth-col"]],
    min_depth = opt[["min-depth"]],
    min_delta = opt[["min-delta"]],
    method = paste0("methylKit_", opt[["methylkit-test"]], "_two-sided"),
    alpha = opt$alpha,
    methylkit_context = opt[["methylkit-context"]],
    candidate_count = nrow(cand),
    testable_count = length(testable),
    depth_filtered_count = sum(!pass),
    tested_count = sum(sapply(results, `[[`, "status") == "ok"),
    failed_count = sum(sapply(results, `[[`, "status") == "failed"),
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
  message(sprintf("Testable sites: %d", length(testable)))
  message(sprintf("Depth filtered: %d", sum(!pass)))
  message(sprintf("Responsive sites: %d (up_in_%s: %d, down_in_%s: %d)",
                  nrow(resp_dt), group2, up_count, group2, down_count))
  message(sprintf("Response output: %s", resp_out))
}

main()
