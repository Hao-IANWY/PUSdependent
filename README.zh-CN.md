# PUS 依赖的假尿苷化 —— 分析代码与数据

<p align="center">
 <a href="README.md">English</a>·<a href="README.zh-CN.md">中文</a>
</p>

本仓库配套以下手稿：

> **Quantitative TRUA-dependent pseudouridine maps across rRNA, mRNA and tRNA in Arabidopsis**
> Shuang Li<sup>1,2†</sup>, Hao Li<sup>1,2†</sup>, Kai He<sup>2,3,4</sup>, Jie Huang<sup>2</sup>, Danmeng Zhu<sup>2,3</sup>, Yu Yu<sup>1,2</sup>, Xuemei Chen<sup>1,2*</sup>. *期刊待定*. 文章链接：待定。

仓库包含论文所使用的分析代码与数据：

- 用于 BID-seq 文库比对和假尿苷（Ψ）位点鉴定的 Snakemake 流程；
- 鉴定 PUS 依赖位点与温度响应位点的下游脚本；
- 每个分析数据集的完整配置；
- 逐图的绘图代码及其输入数据表。

## 仓库结构

| 路径                                                | 内容                                                                                                                                                                                                                    |
| --------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| [`BID-transcript.script/`](BID-transcript.script/) | 比对与位点鉴定流程（`workflows/`）、辅助程序（`bin/`）、下游调用脚本（`downstream/`）。完整文档见自带的 [README](BID-transcript.script/README.md) 与 [CONFIG_GUIDE](BID-transcript.script/CONFIG_GUIDE.md)。        |
| [`configs/`](configs/)                             | 本研究五套分析数据集的完整配置（每套一个 YAML，见下文「数据集与配置」）。                                                                                                                                               |
| [`Figures/`](Figures/)                             | 逐图绘图代码。每个图目录内含 `data_prep.sh`、已提交的输入数据表（`*.tsv`/`*.bed`）、`plot.R` 以及渲染好的 PDF/PNG/SVG，因此不重跑流程即可直接用已提交的数据表重画每一张图。统一样式代码在 `Figures/shared/`。 |

原始 FASTQ 文件与大型中间结果（BAM、mpileup 表）不在本仓库内，见
「数据可用性」。

## 分析流程

BID-seq（bisulfite-induced deletion sequencing，亚硫酸氢盐诱导缺失测
序）文库以 Mock/Treat 成对测序：Treat 文库经亚硫酸氢盐处理，Ψ 与亚硫酸
氢盐形成加合物，逆转录时在 Ψ 位点诱导产生缺失（gap）信号，未修饰的尿苷
则接近零背景；Mock 为未处理对照（方法见何川课题组的 BID-seq：
[Li et al., *Nat Plants* 2025](https://doi.org/10.1038/s41477-024-01894-7)）。

每套数据集由两个 Snakemake 流程之一端到端处理——poly(A)+ RNA 文库用
[`polyA_mapping.smk`](BID-transcript.script/workflows/polyA_mapping.smk)，
tRNA 富集/totalRNA 文库用
[`totalRNA_mapping.smk`](BID-transcript.script/workflows/totalRNA_mapping.smk)（两
个流程都只分析 R2）——依次经过以下阶段：

1. **预处理** —— `trim_galore` 去接头/质量/长度修剪；PCR 去重（poly(A)
   文库用 `fastp --dedup`，totalRNA 文库用 `seqkit rmdup`）；用 `fastp`
   去除测序开头质量较差的碱基（配置中的 `umi_extract`）。
2. **比对** —— 先用 `bowtie2` 比对到第一个配置的参考（poly(A) 文库为
   rRNA 参考序列，tRNA 富集文库为 tRNA/sn/snoRNA 参考序列）；配置了两
   个参考时，未比对上的 reads 再用 `STAR` 比对 Araport11 代表性 cDNA
   序列。
3. **信号提取** —— 把比对得到的 BAM 加工成逐位点计数表：`samFilter`
   链特异性过滤 → `realignGap` 缺失局部重比对 → `samtools calmd` 标注
   错配并过滤错配过多的 reads → `samtools mpileup` 统计每个位点的覆盖
   深度与缺失计数。
4. **gap 校正** —— `adjustGap` 对数据集内全部样品的比对情况联合重新评
   估 gap 信号，按参考分块并行（`adjust_gap.chunk_bp` / `max_sequences`）。
5. **位点鉴定** —— `filter_t_sites` 先按深度阈值筛出可检验的候选 Ψ 位
   点，`bin/Pseudo_callsites.R` 再对各位点做 Mock-vs-Treat 的 Fisher 检
   验，结合逐重复通过标准（配置中的 `filter` 块）判定是否为阳性位点，
   结果写入 `STAT/fisher/<ref_type>/` 下的 `<sample>_sites.Fisher.tsv`
   （全部检验）与 `<sample>_sites.positive.tsv`（阳性位点）。
6. **下游调用** —— 由 positive 位点表鉴定 PUS 依赖位点
   （[`downstream/call_dependent_methylKit.R`](BID-transcript.script/downstream/call_dependent_methylKit.R)）
   与温度响应位点
   （[`downstream/call_temperature_response_methylKit.R`](BID-transcript.script/downstream/call_temperature_response_methylKit.R)）。
   同时，流程还在 `count/<ref_type>/` 下输出逐基因原始表达矩阵（Mock 组
   与全部样本各一份）。

辅助程序 `adjustGap`、`realignGap`、`samFilter`、`cpup` 来自
[y9c/pseudoU-BIDseq](https://github.com/y9c/pseudoU-BIDseq) 项目。

## 位点鉴定统计方法

对每个候选 T（Ψ）位点和每个重复，gap 校正后的 pileup 给出该位点的测序深度（跨越该位点的 reads
数）与 gap 数（在该位点出现缺失的 reads 数），gap 数占测序深度的比例即
为该重复在此位点的 Ψ 比例。这与亚硫酸氢盐测序（BS-seq）对一个胞嘧啶的记录（覆盖度 +
甲基化 reads 数）在形式上完全一致——下文两个差异检验正是因此得以复用
methylKit 的差异甲基化分析框架
（[Akalin et al., *Genome Biol* 2012](https://doi.org/10.1186/gb-2012-13-9-r87)）。

### Ψ 位点鉴定 — [`bin/Pseudo_callsites.R`](BID-transcript.script/bin/Pseudo_callsites.R)

**合并与预过滤。** 对每对 Mock/Treat 样本，各重复的逐位点计数表按位点坐
标（chr、pos、strand）合并，某重复中缺失的位点计零；Mock 与 Treat 任一
组平均深度低于 10 的位点先行剔除。

**组水平检验。** 在两组重复分别求和的 2×2 列联表（gap / 非 gap reads ×
Treat / Mock）上做单侧 Fisher 精确检验（备择方向为 Treat 的缺失频率更
高），p 值经 Benjamini–Hochberg 校正后以 0.01 为阈保留。

**逐重复标准。** 除组水平检验外，每个重复还须独立支持该判定：

- Treat gap 数 ≥ 5；
- Treat gap 比例 ≥ 5%；
- Treat gap 比例不低于配对 Mock 重复的两倍。重复按 Mock 第 *i* 个对
  Treat 第 *i* 个成对比较；两组重复数不等时，每个 Treat 重复改为与
  Mock 组均值比较。

**阳性判定。** 位点在通过组水平检验、且至少 `min_pass_reps` 个重复通过
逐重复标准时判为阳性——三重复数据集为 3/3，两重复的叶绿体数据集为
2/2。阳性表中报告的 ratio 与 depth 仅对通过筛选的重复取平均。

### PUS 依赖位点鉴定 — [`downstream/call_dependent_methylKit.R`](BID-transcript.script/downstream/call_dependent_methylKit.R)

**基本思路。** 若某 Ψ 位点由特定 PUS 酶写入，突变该酶应使修饰消失——
该位点 Treat 文库的 gap 比例应回落至 Mock 背景。脚本因此对每个位点只问
一个问题：*突变体的修饰比例是否显著低于 WT？* 并对该基因的两个独立突
变株系分别作答。

**为什么能用 DNA 甲基化工具来做这个检验。** 每个位点、每个重复上，
BID-seq 只产生两个数：多少条 reads 覆盖该位点、其中多少条带缺失。这与
BS-seq 对一个胞嘧啶位点的记录完全对应：覆盖深度对应测序覆盖度，gap 数
对应甲基化 reads 数。因此把 gap 比例当作“甲基化频率”，这个检验就成了
标准的差异甲基化分析问题，可直接交给 methylKit 完成。

**每个突变株系的判定步骤：**

1. **候选池。** 只考察同一数据集 WT 中鉴定为阳性的位点——在突变体中只
   追问「WT 修饰的位点是否仍被修饰」。
2. **深度门槛。** WT 的各 Treat 重复、突变体的各 Mock 与 Treat 重复在该
   位点的 depth 都必须 ≥ 10，任何一个重复不够就弃掉该位点。这是为了让
   「没测到」无法冒充「修饰消失」（校正表中未覆盖位点记为 depth = 0），
   只保留完整可比的位点。Mock 文库只用于见证覆盖度，真正比较的信号是
   Treat 的 gap 比例。
3. **检验。** 把 WT 与突变体各重复在该位点的覆盖度与 gap 计数交给
   methylKit（`methRead` → `unite` → `calculateDiffMeth`）做组间比较：
   以覆盖度加权，检验 gap 比例是否随基因型（WT 对突变体）改变；模型还
   计入重复之间超出二项预期的额外波动（过离散），避免把正常涨落误判为
   显著。一次要检验成百上千个位点，p 值经 SLIM 做跨位点校正以控制假阳
   性。由于只关心「突变体修饰变少」这一个方向，methylKit 给出的双侧 p
   在方向相符（WT 修饰更高）时折半当单侧使用，方向相反则直接置 1（视
   为不显著）。
4. **效应量。** Δ = WT_Ratio − Mut_Ratio，两项分别是各重复 Treat gap 比
   例的算术平均。该值由输入表直接算得、独立于 methylKit，因此报告的差
   值始终与原始比例一致。
5. **判定。** 单侧 p < 0.05 **且** Δ ≥ 0.05 → 在该株系中依赖。最终只保
   留**两个**株系都依赖的位点，单个株系的特殊背景无法让位点通过；报告
   中突变体比例取两株系均值。

输出：逐株系的 `<allele>.dependent-sites.tsv`、合并的
`<gene>.dependent.tsv`，以及一个汇总 JSON。

### 温度响应位点鉴定 — [`downstream/call_temperature_response_methylKit.R`](BID-transcript.script/downstream/call_temperature_response_methylKit.R)

温度响应位点用同一套 methylKit 流程鉴定，比较对象改为同一 WT 材料在两
个生长温度（27 °C 对 21 °C）下的修饰差异。

**候选池与门槛。** 候选池取两个温度阳性位点的并集（仅在单一温度阳性的
位点正是温度依赖修饰变化的强候选），由调用脚本预先合并。逐重复深度门
槛同样施加于两个温度组的全部 Mock 与 Treat 重复（均须 ≥ 10）——这同
时排除了某温度下转录本不表达的情形，避免其伪装成修饰差异。

**检验与判定。** 与依赖位点不同，温度响应不预设方向，直接使用双侧 p
值：双侧 p < 0.05 且 |Δ| ≥ 0.05 判为温度响应（组比例为各温度 Treat 重
复 gap 比例的算术均值），每个位点另行记录方向（27 °C 下升高或降低）。

输出：`<group2>_vs_<group1>.all-sites.tsv`、`.response.tsv`，以及一个含
上调/下调计数的汇总 JSON。

## 数据集与配置

[`configs/`](configs/) 保存了本研究五套分析数据集的完整 YAML 配置（在工
作树中，每个文件以 `config.yaml` 软链接进各自的流程运行目录）：

| 配置                                                      | 流程             | 数据集                                                                  |
| --------------------------------------------------------- | ---------------- | ----------------------------------------------------------------------- |
| [`polyARNA.yaml`](configs/polyARNA.yaml)                 | polyA_mapping    | poly(A)+ RNA，WT 及七个 TRUA/PUS 基因的 14 个突变体株系，根/叶/花三组织 |
| [`tRNA.yaml`](configs/tRNA.yaml)                         | totalRNA_mapping | tRNA 富集文库，基因型与组织同上                                         |
| [`Temperature.yaml`](configs/Temperature.yaml)           | polyA_mapping    | WT 与两个 trua1b 株系在 21 °C / 27 °C 下的 poly(A)+ RNA               |
| [`Temperature_tRNA.yaml`](configs/Temperature_tRNA.yaml) | totalRNA_mapping | tRNA 富集文库，温度设计同上                                             |
| [`Chloroplast.yaml`](configs/Chloroplast.yaml)           | totalRNA_mapping | 叶绿体 RNA                                                              |

### 配置参数说明

五套配置共用同一套结构：

- **`paths`** —— 辅助程序（`script_dir`、`pseudo_callsites_r`）与比对索
  引（`index_bowtie2` 对应第一个参考；双参考时才需要 `index_star`）的绝
  对路径。
- **`env`** —— 流程脚本所用的 Python 与 R 可执行文件。
- **`pipeline.ref`** —— 一至两个命名参考 FASTA，按比对顺序排列：第一
  个为 bowtie2 阶段，（可选的）第二个为 STAR 阶段。参考名同时作为输出
  目录名（`STAT/<ref_type>/`、`count/<ref_type>/`）。
- **`pipeline.t_site_windows`** —— 候选 T 位点两侧提取的序列上下文碱
  基数。
- **`pipeline.expression`** —— 可选的逐参考原始表达矩阵；`count_mode`
  为 `all_alignments`（比对级计数，适合 rRNA/tRNA 等多比对 RNA）或
  `primary`（read 级计数，适合 mRNA）。
- **`pipeline.adjust_gap`** —— 联合 gap 校正的分块参数：每块
  `chunk_bp`，小转录本按 `max_sequences` 合批。
- **`pipeline.filter`** —— 位点鉴定阈值：
  - `min_depth`：候选池的单重复深度；
  - `min_avg_depth`：组平均深度；
  - `max_padj`：BH 校正 Fisher p 值上限；
  - `min_gap`：每个通过重复的 Treat gap 数；
  - `min_treat_ratio`：每个通过重复的 Treat gap 比例；
  - `ratio_fold`：Treat/Mock 比例倍数；
  - `min_pass_reps`：至少多少个重复通过——三重复数据集为 3，两重复的
    Chloroplast 数据集为 2。
- **`read_layout` / `rmdup` / `trim_galore_length`** —— 双端/单端输
  入、是否去重、修剪后保留的最短 read 长度。
- **`pipeline.umi_extract`** —— 用 `fastp` 从 read 两端做固定长度切除
  （`trim_front`/`trim_tail` 分别为 5′/3′ 端切除的碱基数）。切除长度由
  文库结构决定：本研究的文库在每条 read 上带有固定的 8 nt UMI 与
  `TATAGGG` 序列，同时需去掉 5′ 端测序质量较差的前几个碱基。
  - **poly(A) 文库**：去掉前 6 nt 低质量碱基，故 5′ 端共切除 21 nt
    （8 + 7 + 6）、3′ 端切除 6 nt，双端文库按两端互补设为
    `R1: 6/21, R2: 21/6`。
  - **tRNA 富集/totalRNA 文库**：reads 较短，只去除前 2 个低质量碱基；
    R2 的 5′ 端结构为 `8 bp UMI + TATAGGG + 2 bp`，即 `R2: 17/0`（R1
    镜像为 `0/17`）。
  - 无论 poly(A) 还是 tRNA/totalRNA 文库，比对及下游分析都只使用 R2。
- **`samples`** —— 逐样本的重复定义：
  `samples.<sample>.{mock,treat}.<rep>.{file_index,read1,read2}`，FASTQ
  路径相对于运行目录的 `raw_data/`。每个重复可选 **`stages`** 键限定其
  进入的分支：`stages: [count]` 的重复只参与比对/表达矩阵，不进 Fisher
  位点鉴定分支，因此新增这类重复不会使联合 gap 校正结果失效（默认值
  `[count, fisher]` 两条分支都进）。

### 运行一套数据集

每个运行目录存有各自的 `config.yaml`（软链接自 `configs/`）和存放
FASTQ 的 `raw_data/` 目录，执行：

```bash
snakemake \
 -s /path/to/BID-transcript.script/workflows/polyA_mapping.smk \
 --configfile config.yaml \
 --directory . \
 --cores 200 \
 --rerun-incomplete
```

（totalRNA 数据集改用 `totalRNA_mapping.smk`）。主要输出为
`STAT/mpileup/`、`STAT/adjustGap/`、`STAT/fisher/`（positive 位点表）、
`count/`（表达矩阵）与 `logs/`（逐重复文库统计）。

## 环境依赖

### 比对与位点鉴定流程

- **编排**：[Snakemake](https://snakemake.github.io)。
- **命令行工具**：`trim_galore`（内部调用 cutadapt 与 FastQC）、`fastp`、
  `seqkit`、`bowtie2`、`STAR`、`samtools`、`bedtools`、`csvtk`。
- **Python 3**：[`polars`](https://pola.rs/)、
  [`pysam`](https://pysam.readthedocs.io/)、PyYAML。
- 辅助程序 `adjustGap`、`realignGap`、`samFilter`、`cpup` 已随仓库附带于
  `BID-transcript.script/bin/`，无需另行安装。

### 下游差异检验

- Bioconductor 包
  [`methylKit`](https://bioconductor.org/packages/methylKit/)（PUS 依赖位点
  与温度响应位点鉴定）。

## 数据可用性

本研究生成的原始测序数据已上传至 **[待定]**（链接待定）。每张图所依据
的逐位点处理数据表已随绘图代码一并提交在 [`Figures/`](Figures/) 下。

## 引用

若您使用了本仓库的代码或数据，请引用：

> [引用信息待定]
