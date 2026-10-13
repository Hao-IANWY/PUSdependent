# BID-transcript.script

BID-seq 转录组分析的脚本合集：Snakemake 工作流 + 流程内工具 + 独立下游工具。

## 目录结构

```
BID-transcript.script/
├── README.md                    ← 本文件
├── CONFIG_GUIDE.md              ← config.yaml 编写说明（单/双端、UMI、去重等）
├── config.example.yaml          ← 配置模板，新分析复制后修改
├── workflows/                   ← Snakemake 工作流（snakemake -s 指定）
│   ├── polyA_mapping.smk                polyA 主流程（含分块 adjustGap）
│   ├── totalRNA_mapping.smk             totalRNA 流程
│   └── polyA_mapping-adjustDeletion.smk polyA + deletion 右移校正变体
├── bin/                         ← 流程内工具，即 config 里的 script_dir
│   ├── samFilter                BAM 过滤（C 二进制）
│   ├── cpup                     pileup 计数（C 二进制）
│   ├── realignGap               gap 重比对（Python 脚本，直接执行）
│   ├── adjustGap                gap 校正（Python 脚本，直接执行）
│   ├── adjustDeletion.py        deletion 校正（仅 adjustDeletion 变体使用）
│   ├── combine_mpileup_tsv.py   合并各参考类型的 mpileup tsv
│   ├── count_feature_alignments.py  表达 count 矩阵（polyA/totalRNA 可选）
│   ├── write_adjust_gap_chunks.py   ┐
│   ├── split_mpileup_chunks.awk     ├ 分块 adjustGap 三件套（仅 polyA 主流程）
│   ├── concatenate_adjusted_chunks.py ┘
│   └── Pseudo_callsites.R       Fisher 检验位点调用（config 的 pseudo_callsites_r）
└── downstream/                  ← 独立下游工具，不被任何 smk 调用，手动运行
    ├── call_dependent.py        两基因型比较，调用 PUS 依赖位点
    ├── fit_calibration.py       按校准曲线把 ratio 换算为 fraction
    ├── convert_to_transcript_coord.py  基因组坐标位点表 → 本流程转录组坐标
    └── Calibration_curve.tsv    校准曲线数据（fit_calibration.py 默认输入）
```

## 调用关系

- 工作流内所有工具均通过 config 的 `paths.script_dir` 引用（形如 `{SCRIPT_DIR}/samFilter`），
  **script_dir 必须指向 `bin/` 目录**；R 脚本由独立的 `paths.pseudo_callsites_r` 指定。
- 因此 `bin/` 是整体部署单元：同步到生产拷贝时整个 `bin/` 一起拷，不要拆开。
- 3 个 smk 不硬编码任何工具的相对位置，移动本仓库位置只需改 config 里的路径。

## 运行工作流

```bash
snakemake -s workflows/polyA_mapping.smk --configfile config.yaml --cores 120
```

config.yaml 从 `config.example.yaml` 复制修改，字段说明见 `CONFIG_GUIDE.md`。

## 表达 count 输出

`polyA_mapping.smk` 与 `totalRNA_mapping.smk` 支持可选 `pipeline.expression` 配置。启用后每个 ref_type 固定生成两张原始 count 矩阵：

- `all_alignments`：统计所有 mapped alignment，适用于保守 rRNA/ncRNA，一条 read 可进入多个 feature；
- `primary`：只统计 primary alignment，适用于 mRNA，一条 read 至多进入一个 transcript；
- `count/<ref_type>/Mock.<mode>_counts.tsv`：仅 mock 组样品，沿用既有结构；
- `count/<ref_type>/all_samples.<mode>_counts.tsv`：全部样品（mock 与 treat）。

矩阵包含参考 FASTA 中的全部 feature 和原始整数 count。

两类模式的计数总和语义不同，不能直接相加或直接比较总量。

## 下游工具

### call_dependent.py — PUS 依赖位点调用

以 WT 候选位点 TSV 为基准，比较两个突变株系的 Treat/Mock adjusted TSV：

```bash
python downstream/call_dependent.py \
    --tsv WT_sites.positive.tsv \
    --genotype1 trua1a_1 --genotype2 trua1a_2 \
    --genotype1-treat-files Treat.g1.rep1.adjusted.tsv,Treat.g1.rep2.adjusted.tsv \
    --genotype1-mock-files  Mock.g1.rep1.adjusted.tsv,Mock.g1.rep2.adjusted.tsv \
    --genotype2-treat-files Treat.g2.rep1.adjusted.tsv,Treat.g2.rep2.adjusted.tsv \
    --genotype2-mock-files  Mock.g2.rep1.adjusted.tsv,Mock.g2.rep2.adjusted.tsv \
    --outdir dependent_sites
```

### fit_calibration.py — 校准曲线换算

给 TSV 的 ratio 列按 context 换算出 fraction 列（默认曲线为同目录的
`Calibration_curve.tsv`，可用 `--curve` 覆盖）：

```bash
python downstream/fit_calibration.py --tsv sites.Fisher.tsv --output sites.calibrated.tsv
```

### convert_to_transcript_coord.py — 基因组坐标 → 转录组坐标

本流程把 reads 直接比对到转录组 cDNA 参考，`STAT/` 下所有表的 `chr` 列是**转录本 ID**、
`pos` 是该转录本内的 1-based 坐标。外部流程（在基因组上比对调用）给出的位点表用的是
染色体 + 基因组坐标，两者无法直接比对，需要本脚本换算。

两种模式：

**`--mode tx-col`（默认）** — 按表内自带的转录本 ID 列换算，允许非代表异构体：

```bash
python downstream/convert_to_transcript_coord.py \
    --mode tx-col \
    --input  sites.genome.tsv \
    --gtf    <REPO_ROOT>/Araport11_fixed.gtf \
    --cdna   /home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model \
    --output sites.transcript.tsv
```

追加 4 列：`tx_pos`、`tx_pos_type`（`exon` / `intron_N` / `upstream` / `downstream`）、
`tx_base`、`tx_base_match`。

**`--mode genome`** — **忽略**表内转录本 ID，只用 (chr, pos, strand) 基因组坐标，
在本流程实际使用的**代表转录本**集合（由 `--cdna` 的 `.fai` 定义）里反查：

```bash
python downstream/convert_to_transcript_coord.py \
    --mode genome \
    --input  sites.genome.tsv \
    --gtf    <REPO_ROOT>/Araport11_fixed.gtf \
    --cdna   /home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model \
    --output sites.rep_transcript.tsv
```

追加 8 列：`rep_tx_id`、`rep_tx_pos`、`rep_tx_base`、`rep_base_match`、
`rep_assign`（`unique` / `ambiguous` / `antisense_only` / `none`）、`rep_n_candidates`、
`rep_candidates`（全部候选，`tx:pos` 分号分隔）、`rep_spanning`（外显子反查落空时，
列出 mRNA 跨度覆盖该位置的代表转录本及其类型，用来说明落空原因）。

#### 使用要点

- 只应保留外显子命中（tx-col 模式 `tx_pos_type == exon`；genome 模式 `rep_assign != none`）。
  落在内含子或转录本跨度外的位点在本流程的 cDNA 参考上没有对应位置。
- tx-col 模式下若 `tx_base_match` 为 `FALSE`，说明该行的转录本 ID 与位点实际所在转录本
  不一致（多异构体基因常见），应改用 genome 模式按坐标反查。
- genome 模式同一坐标落在多条代表转录本上时（基因重叠）取主选并标记 `ambiguous`，
  主选规则为优先落在 CDS 内、其次取更长的转录本、最后按 ID 字典序；全部候选见
  `rep_candidates`。
- **坐标语义与 gppy `gtf.py` 完全一致**（已用 23071 个外显子/内含子/上下游探针差分测试，
  0 分歧），因此两者可互相替代。
- 该 GTF 的外显子拼接序列与 `--cdna` 指定的代表 cDNA FASTA 逐碱基一致
  （27526/27562 条；35 条例外为线粒体 ATMG 版本差异与 7 条核基因注释差异）。
