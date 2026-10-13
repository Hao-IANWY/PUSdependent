# bin/ — 流程内工具

`bin/` 是 config 里 `paths.script_dir` 指向的**整体部署单元**：同步到生产拷贝时
整个 `bin/` 一起拷，不要拆开。除 `Pseudo_callsites.R` 由独立的
`paths.pseudo_callsites_r` 指定外，所有工具由 smk 以 `{SCRIPT_DIR}/<工具>` 形式调用。

## 核心链（按数据流顺序）

- `samFilter` — BAM 过滤（C 二进制）
- `cpup` — pileup 计数（C 二进制）
- `realignGap` — gap 重比对（Python 脚本，直接执行）
- `adjustGap` — gap 校正（Python 脚本，直接执行）

## 分块 adjustGap 三件套（仅 polyA 主流程）

- `write_adjust_gap_chunks.py`
- `split_mpileup_chunks.awk`
- `concatenate_adjusted_chunks.py`

## 合并与计数

- `combine_mpileup_tsv.py` — 合并各参考类型的 mpileup tsv
- `count_feature_alignments.py` — 表达 count 矩阵（polyA/totalRNA 可选）

## 文库日志三件套

- `build_library_log.py` — 构建单文库 mapping/deletion 日志 TSV
- `summarize_processing_logs.py` — 合并各文库各阶段读数成长表
- `summarize_deletion_signals.py` — 汇总 cpup gap 缺失信号与含 deletion 的 reads

## 位点调用

- `Pseudo_callsites.R` — Fisher 检验位点调用（config 的 `pseudo_callsites_r`）

## 独立工具

- `adjustDeletion.py` — deletion 校正（当前无 smk 调用）
