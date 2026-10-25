#!/usr/bin/env Rscript
# ============================================================================
# call_dependent_methylKit.R — methylKit 版 dependent 位点判定
#
# 与 call_dependent_stats.R（quasi-binomial / glmmTMB 计数模型版）并列，改用
# methylKit 对修饰差异做统计学检验，做法对齐 Akalin et al. 2012 (Genome Biology
# 13:R87) 的经典工作流：
#
#   methRead()           逐个样本读入甲基化计数
#   unite()              合并样本（complete-case：位点需在所有样本出现）
#   calculateDiffMeth()  差异甲基化检验。methylKit 按样本总数分流：仅当总共
#                        2 个样本（每组 1 个）时 test 被强制改写为 fast.fisher
#                        （Fisher 精确检验）；有生物学重复时（本项目 3v3）走
#                        logReg：以 coverage 为权的二项 logistic 回归，
#                        overdispersion="MN" 估计过离散因子 phi，phi>1 时做
#                        F 检验（deviance/phi ~ F(1, n-2)），phi<=1 退化为卡方。
#
# 分组约定：WT 记作 treatment=1、突变株记作 treatment=0，因此
#     meth.diff = (WT 组加权平均修饰) - (突变组加权平均修饰)，单位百分点。
# 要判定的“依赖位点”（突变株修饰丢失）即 meth.diff > 0 的方向。
#
# 显著性口径（--significance）：
#   q（默认）  用 methylKit calculateDiffMeth(--methylkit-adjust，默认 BH) 的
#              qvalue 判定，同时要求方向为正（meth.diff > 0）。
#   p          旧的原始单侧 p 口径（p_two/2，方向为负时取 1），仅用于复现旧结果。
# qvalue 是双侧 FDR；多重校正默认 BH（2026-10-20 起。此前用 SLIM，其逐线
# pi0 估计在本项目 100~200 位点的候选宇宙上被系统性低估、q 偏松甚至可为
# 负，--methylkit-adjust SLIM 仅保留作复现当日 SLIM 版产物用）。
# 单侧 p 值仍会写入输出表的 p.value 列，供审计和下游参考。
#
# 输出口径严格对齐 call_dependent_stats.R / call_dependent.py：
#   WT_Ratio  = --ratio-col（WT Treat 平均修饰比）
#   Mut_Ratio = 突变株各 Treat 重复 ratio 的简单均值
#   Delta_Ratio(WT-Mut) = WT_Ratio - Mut_Ratio（简单均值差，用于判定）
# methylKit 的 meth.diff 仅用于方向与显著性，不参与 delta 计算。
# 阈值：--alpha 0.05、--min-delta 0.05（fraction）；两个株系都显著才 dependent。
#
# 深度门槛：突变株 Mock 与 Treat 的每一个重复 depth 都必须 >= --min-depth
# （默认 10），WT 每个重复同样要求 >= --min-depth；任一重复不达标即剔除该
# 位点（status="depth_filter"）。adjusted 全量网格对未覆盖位点写 depth=0 且
# ratio=0（防除零约定，非“修饰消失”的有效观测），该门槛保证任何零覆盖或
# 低覆盖重复都不会进入检验，因此无需零覆盖特判逻辑。
# ============================================================================
# 历史文档归档：README.zh-CN.md 原“PUS 依赖位点鉴定”方法描述
# 主文档于 2026-10-25 改用 beta-binomial + LRT；以下保留旧口径供追溯。
# ---------------------------------------------------------------------------
# 基本思路。若某 Ψ 位点由特定 PUS 酶写入，突变该酶应使修饰消失——该位点
# Treat 文库的 gap 比例应回落至 Mock 背景。脚本因此对每个位点只问一个问题：
# 突变体的修饰比例是否显著低于 WT？并对该基因的两个独立突变株系分别作答。
#
# 为什么能用 DNA 甲基化工具来做这个检验。每个位点、每个重复上，BID-seq 只
# 产生两个数：多少条 reads 覆盖该位点、其中多少条带缺失。这与 BS-seq 对一个
# 胞嘧啶位点的记录完全对应：覆盖深度对应测序覆盖度，gap 数对应甲基化 reads
# 数。因此把 gap 比例当作“甲基化频率”，这个检验就成了标准的差异甲基化分析
# 问题，可直接交给 methylKit 完成。
#
# 每个突变株系的判定步骤：
#
#   1. 候选池。只考察同一数据集 WT 中鉴定为阳性的位点——在突变体中只追问
#      “WT 修饰的位点是否仍被修饰”。
#   2. 深度门槛。WT 的各 Treat 重复、突变体的各 Mock 与 Treat 重复在该位点
#      的 depth 都必须 >= 10，任何一个重复不够就弃掉该位点。这是为了让“没测到”
#      无法冒充“修饰消失”（校正表中未覆盖位点记为 depth = 0），只保留完整
#      可比的位点。Mock 文库只用于见证覆盖度，真正比较的信号是 Treat 的 gap
#      比例。
#   3. 检验。把 WT 与突变体各重复在该位点的覆盖度与 gap 计数交给 methylKit
#      （methRead → unite → calculateDiffMeth）做组间比较：以覆盖度加权，检验
#      gap 比例是否随基因型（WT 对突变体）改变；模型还计入重复之间超出二项
#      预期的额外波动（过离散），避免把正常涨落误判为显著。一次要检验成百上千
#      个位点，双侧 p 值经 Benjamini–Hochberg 校正折算成 q 值（双侧 FDR）以
#      控制跨位点假阳性；判定直接使用该 q 值。由于只关心“突变体修饰变少”这
#      一个方向，还要求方向相符（meth.diff > 0，WT 修饰更高），方向相反则一律
#      不通过。
#   4. 效应量。Delta = WT_Ratio − Mut_Ratio，两项分别是各重复 Treat gap 比例
#      的算术平均。该值由输入表直接算得、独立于 methylKit，因此报告的差值始终
#      与原始比例一致。
#   5. 判定。BH q < 0.05、方向为正（meth.diff > 0）且 Delta >= 0.05，则在该
#      株系中依赖。最终只保留两个株系都依赖的位点，单个株系的特殊背景无法让
#      位点通过；报告中突变体比例取两株系均值。
#
# 输出：逐株系的 <allele>.dependent-sites.tsv、合并的 <gene>.dependent.tsv，
# 以及一个汇总 JSON。
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
              help = "Multiple-testing adjustment in calculateDiffMeth: BH (default) or SLIM (legacy; its per-run pi0 estimate is unreliable on small candidate sets)."),
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
# 把 WT / mutant 样本写成 methylKit 泛型格式（chr start end strand coverage freqC）
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
      meth, overdispersion = "MN", adjust = opt[["methylkit-adjust"]],
      test = opt[["methylkit-test"]], mc.cores = opt[["mc-cores"]]
    )
    getData(diff)
  }, error = function(e) {
    message(sprintf("[%s] methylKit failed: %s", name, conditionMessage(e)))
    NULL
  })
}

# ---------------------------------------------------------------------------
# 单株系主流程：读 mock/treat，逐位点算效应量，仅对每重复深度达标的位点跑
# methylKit 检验，返回单侧 p 值。
# ---------------------------------------------------------------------------
call_line_methylkit <- function(name, wt_rows, wt_ratio_mat, wt_depth_mat, mock_files, treat_files) {
  target_keys <- wt_rows$site_id
  allow_dup <- opt[["allow-duplicate-sites"]]

  mock_tables  <- lapply(mock_files,  read_adjusted, target_keys = target_keys, allow_dup = allow_dup)
  treat_tables <- lapply(treat_files, read_adjusted, target_keys = target_keys, allow_dup = allow_dup)

  # complete-case：位点需同时出现在所有 mock 与所有 treat 文件（对齐 call_dependent_stats.R）
  testable <- target_keys
  for (tb in mock_tables)  testable <- intersect(testable, tb$site_id)
  for (tb in treat_tables) testable <- intersect(testable, tb$site_id)

  if (length(testable) == 0) {
    return(list(name = name, results = list()))
  }

  w_idx <- match(testable, wt_rows$site_id)

  # 突变株 Mock/Treat 的 ratio/depth 矩阵。adjusted 是全量位点网格，未覆盖位点
  # depth=0 且 ratio 硬编码为 0，由下面的每重复深度门槛整体挡在检验之外。
  mock_depth_mat  <- do.call(cbind, lapply(mock_tables,  function(tb) tb$depth[match(testable, tb$site_id)]))
  treat_ratio_mat <- do.call(cbind, lapply(treat_tables, function(tb) tb$ratio[match(testable, tb$site_id)]))
  treat_depth_mat <- do.call(cbind, lapply(treat_tables, function(tb) tb$depth[match(testable, tb$site_id)]))

  mut_ratio <- rowMeans(treat_ratio_mat)
  wt_ratio  <- wt_rows$wt_ratio[w_idx]
  delta     <- wt_ratio - mut_ratio

  # 每重复深度门槛：突变株 Mock 与 Treat 的全部重复、WT 的全部重复都须
  # >= min-depth；任何一个重复不达标即剔除（status="depth_filter"）。
  min_d <- opt[["min-depth"]]
  depth_ok <- function(m) (rowSums(m >= min_d) == ncol(m)) & complete.cases(m)
  pass <- depth_ok(mock_depth_mat) & depth_ok(treat_depth_mat) &
          depth_ok(wt_depth_mat[w_idx, , drop = FALSE])

  # 单侧 p 值与 q 值（qvalue 为双侧 FDR，校正方法由 --methylkit-adjust 决定）；
  # 深度不达标的位点留 NA
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
    # q 口径：qvalue < alpha 且方向为正；p 口径保留旧行为。
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
# 输出前缀推断（对齐 call_dependent.py / call_dependent_stats.R）
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
  # 保留 key 列的原始取值，供 write_methyl_sample 写 chr/start/end/strand
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

  # line-level 表（每个株系所有位点）
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
