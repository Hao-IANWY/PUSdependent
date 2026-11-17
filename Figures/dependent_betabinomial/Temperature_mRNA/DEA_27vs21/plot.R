#!/usr/bin/env Rscript
# ============================================================================
# plot.R — DESeq2 差异分析（27°C vs 21°C）× {WT, trua1b_1, trua1b_2} + 火山/MA 图
#
# 铁律：数据只来自 read.table() 读取的 TSV（data_prep.sh 生成），绝不内嵌数据。
#   count_27vs21.tsv : count 矩阵（行=基因，第一列=feature_id；取 Mock rep3/4/5）
#   samples.tsv      : 样品元数据（sample / genotype / temperature）
#
# 可选环境变量：
#   LOGFC_FILTER  |log2FC| 阈值（默认 1）
#   PADJ_FILTER   padj 阈值（默认 0.05）
#   FILTER_TYPE   低表达过滤（默认 "All"：组内 6 个样品 count 全 > 0 才保留）
#   OUTPUT_STEM   组合图文件名主干（默认 "volcano_27vs21"）
#   THEME_R       theme.R 路径（默认脚本目录）
#
# 输出：
#   {OUTPUT_STEM}.pdf/.png/.svg                    三面板组合火山图（WT 最左）
#   {genotype}/volcano_{genotype}_27vs21.*         各基因型独立火山图（offset axes）
#   {genotype}/MAplot_{genotype}_27vs21.*          各基因型 MA plot（闭合边框）
#   {genotype}/{genotype}_27C_vs_21C_DESeq2.csv    各基因型完整结果（统计量+归一化 count）
#   DE_summary.tsv                                 每基因型检验基因数 / Up / Down
# ============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggprism)
  library(ggtext)
})

args_script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
script_dir <- if (length(args_script)) dirname(normalizePath(args_script)) else getwd()
setwd(script_dir)

theme_r      <- Sys.getenv("THEME_R",     unset = "theme.R")
count_tsv    <- Sys.getenv("COUNT_TSV",   unset = "count_27vs21.tsv")
sample_tsv   <- Sys.getenv("SAMPLE_TSV",  unset = "samples.tsv")
logfc_filter <- as.numeric(Sys.getenv("LOGFC_FILTER", unset = "1"))
padj_filter  <- as.numeric(Sys.getenv("PADJ_FILTER",  unset = "0.05"))
filter_type  <- Sys.getenv("FILTER_TYPE", unset = "All")
output_stem  <- Sys.getenv("OUTPUT_STEM", unset = "volcano_27vs21")
csv_stem     <- "27C_vs_21C"    # 文件名避免 °C 字符

if (!file.exists(theme_r)) stop("theme.R not found: ", theme_r)
source(theme_r)

# ---- 画布与字号（base_size ≈ 画布宽 × 2.67；组合图按单面板宽推导）----
FIG_WIDTH   <- 15            # 组合火山图
FIG_HEIGHT  <- 6.2
PANEL_N     <- 3
BASE_COMB   <- round(FIG_WIDTH / PANEL_N * 72 / 27)   # 5 in -> 13 pt
VOL_W <- 6.5; VOL_H <- 6.2                            # 独立火山图
MA_W  <- 6.5; MA_H  <- 6.0                            # MA plot
BASE_FIG    <- round(VOL_W * 72 / 27)                 # 6.5 in -> 17 pt

base_theme <- theme_sciplot_offset(base_size = BASE_COMB)
# element_markdown 不能用 + theme() 合并进 theme_prism 的 element_text（类不同，
# merge 时会被丢弃）；theme_prism 显式定义了 strip.text.x/y，须直接替换槽位，
# 才能把 *trua1b_1* 渲染为 Arial Bold Italic
md_strip <- ggtext::element_markdown(
  family = "Arial", face = "bold", colour = "black",
  size = rel(0.8), margin = margin(6.4, 6.4, 6.4, 6.4))
base_theme$strip.text   <- md_strip
base_theme$strip.text.x <- md_strip
base_theme$strip.text.y <- md_strip
theme_set(base_theme)

# ---- 读取 TSV（唯一数据来源） ----
raw  <- read.table(count_tsv, header = TRUE, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.table(sample_tsv, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
gene_col <- names(raw)[1]

miss <- setdiff(meta$sample, colnames(raw))
if (length(miss) > 0) stop("Samples missing in count TSV: ", paste(miss, collapse = ", "))
message("[sci-plot] Loaded ", nrow(raw), " genes x ", ncol(raw) - 1,
        " samples from ", basename(count_tsv))

TREAT_TEMP <- 27
CTRL_TEMP  <- 21
TREAT_NAME <- "27\u00B0C"
CTRL_NAME  <- "21\u00B0C"
genotypes  <- unique(meta$genotype)          # 面板顺序：WT -> trua1b_1 -> trua1b_2

# ---- 每个基因型独立跑 DESeq2（design ~ condition, 27 vs 21）----
plot_rows <- list()
summ_rows <- list()

for (g in genotypes) {
  sm      <- meta[meta$genotype == g, , drop = FALSE]
  treat_s <- sm$sample[sm$temperature == TREAT_TEMP]
  ctrl_s  <- sm$sample[sm$temperature == CTRL_TEMP]
  if (length(treat_s) < 2 || length(ctrl_s) < 2)
    stop("Genotype ", g, ": need >= 2 replicates per temperature")

  mat <- as.matrix(raw[, c(treat_s, ctrl_s)])
  if (filter_type == "All") {
    keep <- apply(mat, 1, function(x) all(x > 0))
    mat  <- mat[keep, , drop = FALSE]
  } else {
    stop("Unsupported FILTER_TYPE: ", filter_type)
  }
  storage.mode(mat) <- "integer"
  message(sprintf("[sci-plot] %s: filter(All) %d -> %d genes; running DESeq2...",
                  g, nrow(raw), nrow(mat)))

  cond <- factor(rep(c(TREAT_NAME, CTRL_NAME),
                     c(length(treat_s), length(ctrl_s))),
                 levels = c(CTRL_NAME, TREAT_NAME))   # 对照为参考水平
  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData   = data.frame(condition = cond, row.names = c(treat_s, ctrl_s)),
    design    = ~ condition
  )
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("condition", TREAT_NAME, CTRL_NAME))

  resdata <- merge(as.data.frame(res),
                   as.data.frame(counts(dds, normalized = TRUE)),
                   by = "row.names", sort = FALSE)
  names(resdata)[1] <- gene_col

  resdata$significance <- ifelse(
    resdata$log2FoldChange >  logfc_filter & resdata$padj < padj_filter, "Upregulated",
    ifelse(resdata$log2FoldChange < -logfc_filter & resdata$padj < padj_filter,
           "Downregulated", "Not significant"))
  resdata$significance[is.na(resdata$significance)] <- "Not significant"

  dir.create(g, showWarnings = FALSE, recursive = TRUE)
  csv_out <- file.path(g, paste0(g, "_", csv_stem, "_DESeq2.csv"))
  write.table(resdata, csv_out, quote = FALSE, row.names = FALSE, sep = ",")
  message(sprintf("[sci-plot] %s: Up = %d, Down = %d (|log2FC| >= %g, padj <= %g) -> %s",
                  g, sum(resdata$significance == "Upregulated"),
                  sum(resdata$significance == "Downregulated"),
                  logfc_filter, padj_filter, csv_out))

  summ_rows[[g]] <- data.frame(
    genotype     = g,
    genes_tested = sum(!is.na(resdata$padj)),
    up           = sum(resdata$significance == "Upregulated"),
    down         = sum(resdata$significance == "Downregulated"))

  # padj = NA（独立过滤剔除）按 NS 画在 y = 0；padj = 0 的 -log10 封顶 300
  np <- -log10(resdata$padj)
  np[is.na(np)]       <- 0
  np[is.infinite(np)] <- 300
  plot_rows[[g]] <- data.frame(
    genotype       = g,
    log2FoldChange = resdata$log2FoldChange,
    baseMean       = resdata$baseMean,
    negLog10Padj   = np,
    significance   = resdata$significance,
    stringsAsFactors = FALSE)
}

plot_df <- do.call(rbind, plot_rows)
summ_df <- do.call(rbind, summ_rows)
write.table(summ_df, "DE_summary.tsv", quote = FALSE, row.names = FALSE, sep = "\t")

# ---- 面板标签与顺序（用户指定：WT 最左，随后 trua1b_1 / trua1b_2）----
disp_label <- function(g) ifelse(grepl("^[a-z]", g), paste0("*", g, "*"), g)
GENO_LEVELS <- disp_label(genotypes)
plot_df$panel <- factor(disp_label(plot_df$genotype), levels = GENO_LEVELS)

# ---- 轴刻度 helper（刻度必须罩住数据；多处复用） ----
lfc_axis <- function(v) {
  m  <- max(ceiling(max(abs(v), na.rm = TRUE)), 1)
  st <- if (m <= 4) 1 else if (m <= 8) 2 else if (m <= 15) 3 else if (m <= 20) 5 else 10
  li <- ceiling(m / st) * st
  list(limit = li, breaks = seq(-li, li, by = st))
}
np_axis <- function(v) {
  m  <- max(ceiling(max(v, na.rm = TRUE)), 2)
  st <- if (m <= 20) 5 else if (m <= 50) 10 else if (m <= 100) 25 else if (m <= 200) 50 else if (m <= 300) 100 else 200
  li <- ceiling(m / st) * st
  list(limit = li, breaks = seq(0, li, by = st))
}

threshold_y <- -log10(padj_filter)
sig_colors  <- c("Upregulated"     = COLOR_UP,
                 "Downregulated"   = COLOR_DOWN,
                 "Not significant" = COLOR_NS)

# 独立小图的 Up/Down 计数与基因型标识。必须用 annotate()（自带单行 data）：
# 常量参数的 geom_text 若不显式给 data，会继承主数据并按行数重复绘制标签
# （2.2 万个重叠文本框正是此前 geom_label 渲染挂死的根源）
ann_text <- function(x, y, txt, col, hjust, face = "bold") {
  annotate("text", x = x, y = y, label = txt, hjust = hjust, vjust = 1,
           family = "Arial", fontface = face, size = BASE_FIG / .pt,
           color = col)
}

# ============================================================================
# 1) 组合火山图（三面板共享坐标，offset axes）
# ============================================================================
ax <- lfc_axis(plot_df$log2FoldChange); x_limit <- ax$limit; x_breaks <- ax$breaks
ay <- np_axis(plot_df$negLog10Padj);    y_limit <- ay$limit; y_breaks <- ay$breaks

ann_counts <- do.call(rbind, lapply(genotypes, function(g) {
  d <- plot_df[plot_df$genotype == g, ]
  data.frame(
    panel = disp_label(g),
    x     = c(-x_limit, x_limit),
    y     = c(y_limit, y_limit),
    side  = c("down", "up"),
    label = c(sprintf("Down: %d", sum(d$significance == "Downregulated")),
              sprintf("Up: %d",   sum(d$significance == "Upregulated"))))
}))
# 关键：注记图层的 panel 必须与主数据同为 factor 且 levels 一致；
# 字符型 panel 会把 facet 变量降级为字符排序（'*' < 字母，WT 会被排到最后）
ann_counts$panel <- factor(ann_counts$panel, levels = GENO_LEVELS)
ann_down <- ann_counts[ann_counts$side == "down", ]
ann_up   <- ann_counts[ann_counts$side == "up",   ]

p <- ggplot(plot_df, aes(x = log2FoldChange, y = negLog10Padj, color = significance)) +
  geom_point(size = 1.5, alpha = 0.85) +
  geom_hline(yintercept = threshold_y, linetype = "dashed",
             linewidth = 0.6, color = "black") +
  geom_vline(xintercept = c(-logfc_filter, logfc_filter), linetype = "dashed",
             linewidth = 0.6, color = "black") +
  geom_label(data = ann_down, aes(x = x, y = y, label = label),
             hjust = 0, vjust = 1, family = "Arial", fontface = "bold",
             size = BASE_COMB / .pt, color = "#3C5488", fill = "white",
             linewidth = 0, inherit.aes = FALSE, show.legend = FALSE) +
  geom_label(data = ann_up, aes(x = x, y = y, label = label),
             hjust = 1, vjust = 1, family = "Arial", fontface = "bold",
             size = BASE_COMB / .pt, color = "#E64B35", fill = "white",
             linewidth = 0, inherit.aes = FALSE, show.legend = FALSE) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(name = NULL,
                     breaks = c("Upregulated", "Downregulated", "Not significant"),
                     values = sig_colors) +
  labs(x = sprintf("Log2 Fold Change (%s / %s)", TREAT_NAME, CTRL_NAME),
       y = "-Log10 Adjusted p-value") +
  scale_x_continuous(limits = c(-x_limit, x_limit), breaks = x_breaks,
                     expand = expansion(mult = c(0.03, 0.06)),
                     guide = ggprism::guide_prism_offset()) +
  scale_y_continuous(limits = c(0, y_limit), breaks = y_breaks,
                     expand = expansion(mult = c(0.03, 0.06)),
                     guide = ggprism::guide_prism_offset()) +
  theme(
    legend.position = "top",
    legend.title    = element_text(family = "Arial", face = "bold", colour = "black"),
    legend.text     = element_text(family = "Arial", face = "bold", colour = "black"),
    plot.margin     = margin(t = 12, r = 14, b = 12, l = 12, unit = "pt")
  ) +
  coord_cartesian(clip = "off")

save_sciplot(p, output_stem, width = FIG_WIDTH, height = FIG_HEIGHT)
message("[sci-plot] Wrote combined volcano: ", output_stem, ".pdf/.png/.svg")

# ============================================================================
# 2) 每基因型子文件夹：独立火山图 + MA plot
# ============================================================================
# 切到小画布字号的全局主题（组合图已保存，theme_set 不再影响它）。
# 注意：不能对 plot 用 + theme_sciplot_offset() —— 全局 strip 槽位是
# element_markdown，与 theme_prism 的 element_text 在 merge 时直接类冲突报错
theme_set(theme_sciplot_offset(base_size = BASE_FIG))
for (g in genotypes) {
  d <- plot_df[plot_df$genotype == g, ]
  up_n   <- sum(d$significance == "Upregulated")
  down_n <- sum(d$significance == "Downregulated")
  g_face <- if (grepl("^[a-z]", g)) "bold.italic" else "bold"   # 基因名斜体

  # ---- 2a) 独立火山图（offset axes，自适应轴） ----
  # 注记 y 位置统一向面板内侧缩 4%，避免文本贴着刻度上限
  ax <- lfc_axis(d$log2FoldChange)
  ay <- np_axis(d$negLog10Padj)
  y_inset <- ay$limit - 0.04 * ay$limit
  p_vol <- ggplot(d, aes(x = log2FoldChange, y = negLog10Padj, color = significance)) +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_hline(yintercept = threshold_y, linetype = "dashed",
               linewidth = 0.6, color = "black") +
    geom_vline(xintercept = c(-logfc_filter, logfc_filter), linetype = "dashed",
               linewidth = 0.6, color = "black") +
    ann_text(-ax$limit, y_inset, sprintf("Down: %d", down_n), "#3C5488", 0) +
    ann_text( ax$limit, y_inset, sprintf("Up: %d",   up_n),   "#E64B35", 1) +
    ann_text( 0,        y_inset, g, "black", 0.5, face = g_face) +
    scale_color_manual(name = NULL,
                       breaks = c("Upregulated", "Downregulated", "Not significant"),
                       values = sig_colors) +
    labs(x = sprintf("Log2 Fold Change (%s / %s)", TREAT_NAME, CTRL_NAME),
         y = "-Log10 Adjusted p-value") +
    scale_x_continuous(limits = c(-ax$limit, ax$limit), breaks = ax$breaks,
                       expand = expansion(mult = c(0.03, 0.12)),
                       guide = ggprism::guide_prism_offset()) +
    scale_y_continuous(limits = c(0, ay$limit), breaks = ay$breaks,
                       expand = expansion(mult = c(0.03, 0.06)),
                       guide = ggprism::guide_prism_offset()) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(
      legend.position = "top",
      legend.title    = element_text(family = "Arial", face = "bold", colour = "black"),
      legend.text     = element_text(family = "Arial", face = "bold", colour = "black"),
      plot.margin     = margin(t = 12, r = 14, b = 12, l = 12, unit = "pt")
    ) +
    coord_cartesian(clip = "off")
  save_sciplot(p_vol, file.path(g, paste0("volcano_", g, "_27vs21")),
               width = VOL_W, height = VOL_H)

  # ---- 2b) MA plot（闭合边框主题） ----
  m <- d[is.finite(d$baseMean) & d$baseMean > 0, ]
  m$log10BaseMean <- log10(m$baseMean)
  ax  <- lfc_axis(m$log2FoldChange)
  xmx <- ceiling(max(m$log10BaseMean))
  y_inset <- ax$limit - 0.04 * 2 * ax$limit
  p_ma <- ggplot(m, aes(x = log10BaseMean, y = log2FoldChange, color = significance)) +
    geom_point(size = 1.8, alpha = 0.85) +
    geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.6, color = "black") +
    geom_hline(yintercept = c(-logfc_filter, logfc_filter), linetype = "dashed",
               linewidth = 0.6, color = "black") +
    ann_text(0,          y_inset, sprintf("Down: %d", down_n), "#3C5488", 0) +
    ann_text(xmx,        y_inset, sprintf("Up: %d",   up_n),   "#E64B35", 1) +
    ann_text(xmx / 2,    y_inset, g, "black", 0.5, face = g_face) +
    scale_color_manual(name = NULL,
                       breaks = c("Upregulated", "Downregulated", "Not significant"),
                       values = sig_colors) +
    labs(x = "Log10 (Base Mean)",
         y = sprintf("Log2 Fold Change (%s / %s)", TREAT_NAME, CTRL_NAME)) +
    scale_x_continuous(limits = c(0, xmx), breaks = 0:xmx,
                       expand = expansion(mult = c(0.02, 0.10))) +
    scale_y_continuous(limits = c(-ax$limit, ax$limit), breaks = ax$breaks,
                       expand = expansion(mult = c(0.03, 0.06))) +
    theme_sciplot(base_size = BASE_FIG) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
    theme(
      legend.position = "top",
      legend.title    = element_text(family = "Arial", face = "bold", colour = "black"),
      legend.text     = element_text(family = "Arial", face = "bold", colour = "black"),
      plot.margin     = margin(t = 12, r = 14, b = 12, l = 12, unit = "pt")
    ) +
    coord_cartesian(clip = "off")
  save_sciplot(p_ma, file.path(g, paste0("MAplot_", g, "_27vs21")),
               width = MA_W, height = MA_H)

  message(sprintf("[sci-plot] %s: volcano_%s_27vs21.* + MAplot_%s_27vs21.* (Up %d / Down %d)",
                  g, g, g, up_n, down_n))
}

message("[sci-plot] Done: combined volcano + per-genotype volcano/MAplot + CSVs + DE_summary.tsv")
