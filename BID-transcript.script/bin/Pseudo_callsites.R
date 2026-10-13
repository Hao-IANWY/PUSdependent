# 输入 cpup 格式的矩阵文件，对输入的文件进行Fisher检验
# chr	pos	strand	context	ratio	depth	gap

# 默认情况下先进行一轮初筛，Mock/Treat 组位点平均深度需要大于等于10
# 要求在组间Fisher需要显著（多个生物学重复的总和显著）
# 对每一个处理/对照组：
# Treat组上缺失信号要大于等于5
# Treat组上缺失信号比例默认要大于等于5%，可通过 --min-treat-ratio 调整
# Treat组上缺失信号比例是Mock组的两倍

library(dplyr)
library(parallel)
library(optparse)

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input files (e.g., Mock:Mock,Treat:Treat) Mock要在Treat前", dest = "input"),
  make_option(c("-n", "--name"), type = "character", help = "Sample names (e.g., Mock:Mock,Treat:Treat) Mock要在Treat前", dest = "name"),
  make_option(c("-k", "--keepall"), action = "store_true", default = FALSE,
              help = "是否保留所有位点信息 (默认不保留)", dest = "keepall"),
  make_option(c("-r", "--min-treat-ratio"), type = "double", default = 0.05,
              help = "每个Treat重复通过筛选所需的最小ratio阈值 (默认: 0.05)", dest = "min_treat_ratio"),
  make_option(c("--min-avg-depth"), type = "double", default = 10,
              help = "Mock/Treat 组平均 depth 阈值 (默认: 10)", dest = "min_avg_depth"),
  make_option(c("--max-padj"), type = "double", default = 0.01,
              help = "BH 校正后 p 值阈值 (默认: 0.01)", dest = "max_padj"),
  make_option(c("--min-gap"), type = "double", default = 5,
              help = "每个Treat重复通过筛选所需的最小gap阈值 (默认: 5)", dest = "min_gap"),
  make_option(c("--ratio-fold"), type = "double", default = 2,
              help = "Treat ratio 需大于等于 Mock ratio 的倍数 (默认: 2)", dest = "ratio_fold"),
  make_option(c("--min-pass-reps"), type = "integer", default = 2,
              help = "位点判定为阳性所需的最少通过重复数 (默认: 2)", dest = "min_pass_reps"),
  make_option(c("-c", "--read-cores"), type = "integer", default = max(1, detectCores() - 1),
              help = "读取输入文件时使用的并行进程数 (默认: detectCores() - 1)", dest = "read_cores"),
  make_option(c("-o", "--output"), type = "character", help = "Output file path", dest = "output")
)

# 解析参数
opt <- parse_args(OptionParser(option_list = option_list))

# 检查参数
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

message("========== 参数设置 ==========")
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
message("==============================")

# 注意Mock组需要在Treat组之前
# Input <- "Mock.WT.rep1_AvgMod.tsv:Mock.WT.rep2_AvgMod.tsv,Treat.WT.rep1_AvgMod.tsv:Treat.WT.rep2_AvgMod.tsv"
file_group <- strsplit(Input, ",")[[1]]
files <- unlist(lapply(file_group, function(x) strsplit(x, ":")[[1]]))

# Input_names <- "Mock.WT.rep1:Mock.WT.rep2,Treat.WT.rep1:Treat.WT.rep2"  
names_group <- strsplit(Input_names, ",")[[1]]
sample_names <- unlist(lapply(names_group, function(x) strsplit(x, ":")[[1]]))
# 提取样本对列表
mock_sample_names <- strsplit(names_group[1], ":")[[1]]
treat_sample_names <- strsplit(names_group[2], ":")[[1]]

# 定义读取函数
read_file <- function(file, sample_name) {
  # 检查文件是否存在
  if (!file.exists(file)) {
    stop(sprintf("❌ 读取失败：找不到文件 %s（样本：%s）", file, sample_name))
  }
  data <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data <- data[, c("chr", "pos", "strand", "context", "ratio", "depth", "gap")]
  # 强制数值列类型
  numeric_cols <- c("ratio", "depth", "gap")
  data[numeric_cols] <- lapply(data[numeric_cols], as.numeric)
  # 重命名列名为样本相关
  colnames(data)[5:7] <- paste0(sample_name, c("_ratio", "_depth", "_gap"))
  return(data)
}
# 读取所有文件
# data_list <- mapply(read_file, files, sample_names, SIMPLIFY = FALSE)
data_list <- mclapply(seq_along(files), function(i) 
  read_file(files[i], sample_names[i]), mc.cores = read_cores)
# 根据矩阵前四列进行聚合
merged_data <- Reduce(function(x, y) merge(x, y, by = c("chr", "pos", "strand", "context"), all = TRUE), data_list)
# 空缺值补为 0 (仅替换数值型列的NA)
numeric_cols <- sapply(merged_data, is.numeric)
merged_data[numeric_cols][is.na(merged_data[numeric_cols])] <- 0

# 根绝Treat/Mock分组统计平均值
Mock_depth <- paste0(strsplit(names_group[1], ":")[[1]], "_depth")
Treat_depth <- paste0(strsplit(names_group[2], ":")[[1]], "_depth")
Mock_ratio <- paste0(strsplit(names_group[1], ":")[[1]], "_ratio")
Treat_ratio <- paste0(strsplit(names_group[2], ":")[[1]], "_ratio")
Mock_gap <- paste0(strsplit(names_group[1], ":")[[1]], "_gap")
Treat_gap <- paste0(strsplit(names_group[2], ":")[[1]], "_gap")
# 计算mock组 平均depth
merged_data$Mock_avg_depth <- rowMeans(as.data.frame(merged_data[, Mock_depth]))
# 计算treat组 平均depth
merged_data$Treat_avg_depth <- rowMeans(as.data.frame(merged_data[, Treat_depth]))
# 计算mock组 平均ratio
merged_data$Mock_avg_ratio <- rowMeans(as.data.frame(merged_data[, Mock_ratio]))
# 计算treat组 平均ratio
merged_data$Treat_avg_ratio <- rowMeans(as.data.frame(merged_data[, Treat_ratio]))
# 计算treat组 平均gap
merged_data$Treat_avg_gap <- rowMeans(as.data.frame(merged_data[, Treat_gap]))

# keepall为FALSE时，执行筛选流程，否则保留所有位点信息
if (!keepall) {
	# 去掉平均depth小于 min_avg_depth 的位点
	merged_data <- merged_data[merged_data$Mock_avg_depth >= min_avg_depth & merged_data$Treat_avg_depth >= min_avg_depth, ]
}

# One-sided Fisher’s Exact Test
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

# 每一个重复生成独立的Pass列（检验每一个位点是否在一个重复上通过筛选）
# Mock/Treat 重复数一致时，保留原始的 Mock[i] vs Treat[i] 成对比较逻辑；
# 重复数不一致时，用 Mock 组平均信号作为背景，分别和每个 Treat 重复比较。
if (length(mock_sample_names) == length(treat_sample_names)) {
  pass_treat_ratio <- Treat_ratio
  pass_treat_depth <- Treat_depth

  for (i in seq_along(mock_sample_names)) {
    mock_name <- mock_sample_names[i]
    treat_name <- treat_sample_names[i]
    pass_col <- paste0("Rep", i, "_pass")
    pvalue_col <- paste0("Rep", i, "_pvalue")
    padj_col <- paste0("Rep", i, "_p.adj")

    # 执行 Fisher 检验并保存 p 值
    mock_depth <- merged_data[paste0(mock_name, "_depth")]
    mock_gap <- merged_data[paste0(mock_name, "_gap")]
    treat_depth <- merged_data[paste0(treat_name, "_depth")]
    treat_gap <- merged_data[paste0(treat_name, "_gap")]
    merged_data[[pvalue_col]] <- apply(cbind(mock_depth - mock_gap, mock_gap,
                                              treat_depth - treat_gap, treat_gap), 1, get_p)
    merged_data[[padj_col]] <- p.adjust(merged_data[[pvalue_col]], method = "BH")

    merged_data[[pass_col]] <- ifelse(
      as.numeric(merged_data[[paste0(treat_name, "_gap")]]) >= min_gap & # Treat组上缺失信号要大于等于 min_gap
      as.numeric(merged_data[[paste0(treat_name, "_ratio")]]) >= min_treat_ratio & # Treat组上缺失信号比例阈值
      # as.numeric(merged_data[[padj_col]]) < max_padj  &  # 组内显著性判断(是否要求在一个重复内，Treat比Mock更显著)
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

pass_columns <- grep("_pass$", colnames(merged_data), value = TRUE)  # 提取所有Pass列
# 生成总Pass列：所有重复的pass列均为1时返回1，否则0
# merged_data$Pass <- as.integer(rowSums(merged_data[pass_columns]) == ncol(merged_data[pass_columns]))
# 生成总Pass列：所有重复的pass列总和 >= min_pass_reps 时返回 1
merged_data$Pass <- as.integer(rowSums(merged_data[pass_columns]) >= min_pass_reps)

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
