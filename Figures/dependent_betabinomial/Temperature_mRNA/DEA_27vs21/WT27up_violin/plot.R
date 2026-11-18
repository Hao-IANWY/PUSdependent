#!/usr/bin/env Rscript
# ============================================================================
# plot.R — WT 27°C 上调基因在 {WT, trua1b_1, trua1b_2} 三组 27°C 表达量的
#          提琴图 + 显著性分析（仅 WT vs trua1b_1、WT vs trua1b_2）
#
# 铁律：数据只来自 read.table() 读取的 TSV（data_prep.sh 生成），绝不内嵌数据。
#   WT27up_expr_27C.tsv : 长表（genotype / sample / feature_id / tpm 等）
#                         tpm = transcripts per million
#
# 展示：WT 中 27°C 显著上调的 891 个基因（padj<0.05, log2FC>=1），
#       在三个基因型 27°C 条件下的表达量分布（log2[TPM + 1]）
#
# 显著性分析设计（easyplot_analysis.R 受保护助手）：
#   单元 unit   = 基因（feature_id）；配对值 = 该基因 3 个重复的 log2(TPM+1) 均值
#   设计 design = paired（同一批 891 个基因跨基因型天然配对）
#   检验        = 配对 t 检验（双侧；n=891 对，t 近似稳健），
#                 效应量 = 均值差（突变体 - WT）+ 95% CI
#   比较族      = {WT vs trua1b_1, WT vs trua1b_2}，Holm 校正（FWER）
#                 （按用户要求，两个突变体之间不做比较）
#
# 可选环境变量：
#   DATA_TSV     数据文件（默认 "WT27up_expr_27C.tsv"）
#   THEME_R      theme.R 路径（默认上级目录 ../theme.R）
#   ANALYSIS_R   easyplot_analysis.R 路径（默认 sci-plot skill 目录）
#
# 输出：WT27up_expression_violin.pdf/.png/.svg
#       WT27up_violin_stats.tsv（配对检验完整统计结果）
# ============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggprism)
  library(ggtext)
})

args_script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(args_script)) dirname(normalizePath(args_script)) else getwd()
setwd(script_dir)

theme_r    <- Sys.getenv("THEME_R",    unset = file.path("..", "theme.R"))
data_tsv   <- Sys.getenv("DATA_TSV",   unset = "WT27up_expr_27C.tsv")
analysis_r <- Sys.getenv("ANALYSIS_R",
                         unset = "/home/HL/.agents/skills/sci-plot/scripts/easyplot_analysis.R")

if (!file.exists(theme_r)) stop("theme.R not found: ", theme_r)
source(theme_r)
if (!file.exists(analysis_r)) stop("easyplot_analysis.R not found: ", analysis_r)
source(analysis_r)

# ---- 画布与字号（base_size ≈ 画布宽 × 2.67：6 in -> 16 pt）----
FIG_WIDTH  <- 6
FIG_HEIGHT <- 5.5
BASE_SIZE  <- round(FIG_WIDTH * 72 / 27)

# element_markdown 不能用 + theme() 合并进 theme_prism 的 element_text（类不同，
# merge 时会被丢弃）；axis.text.x 与 volcano plot.R 的 strip.text 一样直接替换槽位，
# 才能把 *trua1b_1* / *trua1b_2* 渲染成 Arial Bold Italic
base_theme <- theme_sciplot(base_size = BASE_SIZE)
base_theme$axis.text.x <- ggtext::element_markdown(
  family = "Arial", face = "bold", colour = "black",
  margin = ggplot2::margin(t = 12, unit = "pt"))
# y_axis_cover 把面板顶钉在最高刻度，默认 5.5pt 上边距会令顶边框贴边
# （实测 0.04 in < 0.1 in），加宽四边
base_theme$plot.margin <- ggplot2::margin(t = 16, r = 10, b = 10, l = 10,
                                          unit = "pt")
theme_set(base_theme)

# ---- 读取 TSV（唯一数据来源）----
df <- read.table(data_tsv, header = TRUE, sep = "\t",
                 stringsAsFactors = FALSE, check.names = FALSE)

genotypes <- c("WT", "trua1b_1", "trua1b_2")
df[["genotype"]] <- factor(df[["genotype"]], levels = genotypes)
if (any(is.na(df[["genotype"]]))) {
  stop("Unexpected genotype in data: ",
       paste(setdiff(unique(df[["genotype"]]), genotypes), collapse = ", "))
}
stopifnot(is.numeric(df[["tpm"]]))

# 表达量：log2(TPM + 1)
df[["expression"]] <- log2(df[["tpm"]] + 1)
y_data <- max(df[["expression"]])

# ---- 显著性分析：基因配对 t 检验（WT vs 各突变体），Holm 校正 ----
# 每基因 3 重复取均值作为该基因的配对测量值
gene_means <- aggregate(expression ~ feature_id + genotype, data = df, FUN = mean)
wide <- reshape(gene_means, idvar = "feature_id", timevar = "genotype",
                direction = "wide")
names(wide) <- sub("^expression[.]", "", names(wide))

COMPARISONS <- c("trua1b_1", "trua1b_2")   # 仅 WT vs 各突变体
analysis <- lapply(COMPARISONS, function(g) {
  d <- data.frame(
    unit       = rep(wide[["feature_id"]], 2L),
    genotype   = rep(c("WT", g), each = nrow(wide)),
    expression = c(wide[["WT"]], wide[[g]])
  )
  easyplot_analyze_two_groups(
    data   = d, group = "genotype", value = "expression", unit = "unit",
    order  = c("WT", g),            # 效应量 = 突变体 - WT
    design = "paired", missing = "error", confidence = 0.95
  )
})

adj <- easyplot_adjust_pvalues(
  pvalues = vapply(analysis, function(r) r$test$p_value, numeric(1)),
  labels  = sprintf("WT vs %s", COMPARISONS),
  method  = "holm",
  family  = "Expression of WT 27C-up genes: WT vs each trua1b mutant (paired by gene)"
)

# ---- 统计结果 TSV ----
stats_out <- do.call(rbind, lapply(seq_along(analysis), function(i) {
  r <- analysis[[i]]; eff <- r$effect; tst <- r$test; a <- adj$results[[i]]
  data.frame(
    comparison      = a$comparison,
    test_method     = tst$method,
    n_pairs         = r$n_pairs,
    mean_diff_log2_tpm = signif(eff$estimate, 4),
    ci95_low_log2_tpm   = signif(eff$ci_low, 4),
    ci95_high_log2_tpm  = signif(eff$ci_high, 4),
    t_statistic     = signif(tst$statistic, 4),
    df              = round(tst$df, 1),
    p_value         = format(tst$p_value, digits = 3, scientific = TRUE),
    p_holm          = format(a$p_adjusted, digits = 3, scientific = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.table(stats_out, "WT27up_violin_stats.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
message("[sci-plot] Wrote: WT27up_violin_stats.tsv")

# 括号使用用户指定的 ** 标注；按 Holm 校正后 p < 0.01 判定，
# 精确原始/校正 p 值保存在 WT27up_violin_stats.tsv 中。
fmt_p <- function(p) {
  if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

# 配色：WT 用深红色突出；两个突变体使用相近的浅珊瑚色，保持同组视觉关系。
GENO_COLORS <- c(
  WT       = "#C43C39CC",
  trua1b_1 = "#F2A6A6CC",
  trua1b_2 = "#F5B9ACCC"
)

# ---- 显著性括号几何 ----
# 括号与 ** 也参与 y 轴范围计算，避免注释层把面板顶部额外撑出空白。
stub     <- 0.45                                   # 括号两端竖线的下探深度
h1       <- y_data + 0.8                           # WT vs trua1b_1 括号高度
h2       <- y_data + 2.1                           # WT vs trua1b_2（更高一档）
brk_df <- rbind(
  data.frame(cmp = "WT vs trua1b_1", x = c(1, 1, 2, 2),
             y = c(h1, h1 + stub, h1 + stub, h1)),
  data.frame(cmp = "WT vs trua1b_2", x = c(1, 1, 3, 3),
             y = c(h2, h2 + stub, h2 + stub, h2))
)
lab_df <- data.frame(
  cmp   = c("WT vs trua1b_1", "WT vs trua1b_2"),
  x     = c(1.5, 2.0),
  y     = c(h1, h2) + stub + 0.30,
  label = vapply(seq_along(analysis),
                 function(i) fmt_p(adj$results[[i]]$p_adjusted), character(1))
)

# 用所有实际绘制对象计算刻度：最高刻度覆盖括号与文字，同时不留下多余顶部空间。
ax <- y_axis_cover(c(df[["expression"]], brk_df[["y"]], lab_df[["y"]]))

p <- ggplot(df, aes(x = .data[["genotype"]], y = .data[["expression"]],
                    fill = .data[["genotype"]])) +
  geom_violin(scale = "width", width = 0.82, adjust = 1.05,
              linewidth = 0.6, colour = "#333333", show.legend = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, linewidth = 0.5,
               fill = "white", colour = "#333333", show.legend = FALSE) +
  geom_path(data = brk_df,
            aes(x = x, y = y, group = cmp),
            inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = lab_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            family = "Arial", fontface = "bold", size = BASE_SIZE / 3.0) +
  labs(
    x = NULL,
    y = "WT 27\u00B0C-up genes\nlog2(TPM + 1)"
  ) +
  scale_fill_manual(values = GENO_COLORS) +
  scale_x_discrete(labels = c("WT", "*trua1b_1*", "*trua1b_2*")) +
  scale_y_continuous(breaks = ax$breaks, expand = ax$expand)

save_sciplot(p, "WT27up_expression_violin", width = FIG_WIDTH, height = FIG_HEIGHT)
message("[sci-plot] Done: WT27up_expression_violin.{pdf,png,svg}")
