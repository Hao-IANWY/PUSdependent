# PUS 依赖的假尿苷化 —— 分析代码与数据

<p align="center">
  <a href="README.md">English</a> · <a href="README.zh-CN.md">中文</a>
</p>

本仓库包含了以下手稿所使用的分析代码与处理数据：

> **[Title TBD]**
> [Authors TBD]. *Journal TBD*. 文章链接：待定。

仓库主要内容包括：

- **Snakemake 流程**：用于 BID-seq 文库比对与假尿苷（Ψ）位点鉴定的自动化分析流程。
- **下游分析脚本**：用于鉴定 PUS 依赖性位点与温度响应性 Ψ 位点的统计方法。
- **分析配置**：研究中所有数据集的 YAML 格式配置文件。
- **可视化代码**：用于绘图的脚本与输入数据表，确保所有图表均可复现。

## 仓库结构

| 路径 | 说明 |
|---|---|
| [`BID-transcript.script/`](BID-transcript.script/) | 核心分析流程 (`workflows/`)、辅助程序 (`bin/`) 及下游统计分析脚本 (`downstream/`)。详见自带的 [README](BID-transcript.script/README.md) 与 [CONFIG_GUIDE](BID-transcript.script/CONFIG_GUIDE.md)。 |
| [`configs/`](configs/) | 各数据集的配置文件（每套数据集一个 YAML）。详见下文「数据集与配置」。 |
| [`Figures/`](Figures/) | 绘图脚本与原始数据。每个子目录内含 `data_prep.sh`、所需的输入表 (`*.tsv`/`*.bed`)、`plot.R` 脚本及渲染输出 (PDF/PNG/SVG)。 |

*注：原始 FASTQ 文件及大型中间结果（如 BAM、mpileup 表）未包含在仓库中；详见「数据可用性」。*

## 分析流程

BID-seq（亚硫酸氢盐诱导缺失测序）文库成对测序（Mock/Treat）。在 Treat 文库中，亚硫酸氢盐处理将 Ψ 转化为 Ψ-亚硫酸氢盐加合物，逆转录时会在 Ψ 位点诱导产生特异性缺失，而未处理的 Mock 文库背景极低（[Li et al., *Nat Plants* 2025](https://doi.org/10.1038/s41477-024-01894-7)）。

数据集采用以下两种 Snakemake 流程之一进行端到端处理：
- [`polyA_mapping.smk`](BID-transcript.script/workflows/polyA_mapping.smk)：用于 poly(A)+ RNA 文库。
- [`totalRNA_mapping.smk`](BID-transcript.script/workflows/totalRNA_mapping.smk)：用于 tRNA 富集文库/总 RNA 文库。

流程均只分析 R2 reads，主要阶段如下：

1. **预处理**：使用 `trim_galore` 进行接头去除与质量过滤，使用 `fastp --dedup`（针对 poly(A)）或 `seqkit rmdup`（针对总 RNA）进行 PCR 去重。通过 `fastp` 对 read 末端进行 UMI/接头区域的固定裁剪，参数在配置文件中根据末端定义 (`umi_extract`)。
2. **比对**：首先使用 `bowtie2` 比对到一级参考序列（poly(A) 文库为 rRNA；tRNA 富集文库为 tRNA/sn/snoRNA）。若配置了二级参考序列，未比对上的 reads 将进一步使用 `STAR` 比对到 Araport11 代表性 cDNA 模型。
3. **信号提取**：将 BAM 文件处理为位点特异性计数表：`samFilter`（链特异性过滤）→ `realignGap`（缺失局部重比对）→ `samtools calmd`（错配标注与过滤）→ `samtools mpileup`（位点深度与缺失计数）。
4. **缺失（Gap）校正**：`adjustGap` 在所有重复样本间联合重新评估缺失信号，并按参考序列块并行处理 (`adjust_gap.chunk_bp` / `max_sequences`)。
5. **位点鉴定**：首先按深度过滤候选 Ψ 位点 (`filter_t_sites`)。随后，`bin/Pseudo_callsites.R` 根据配置文件中的标准执行重复内（replicate-aware）的 Mock-vs-Treat Fisher 检验，结果输出至 `STAT/fisher/<ref_type>/` (例如 `<sample>_sites.Fisher.tsv` 与 `<sample>_sites.positive.tsv`)。
6. **下游分析**：分别使用 [`downstream/call_dependent_betabinomial.R`](BID-transcript.script/downstream/call_dependent_betabinomial.R) 与 [`downstream/call_temperature_response_betabinomial.R`](BID-transcript.script/downstream/call_temperature_response_betabinomial.R) 鉴定 PUS 依赖位点与温度响应位点。流程同时在 `count/<ref_type>/` 下生成逐基因原始表达矩阵。

辅助程序 `adjustGap`、`realignGap`、`samFilter` 与 `cpup` 来自 [y9c/pseudoU-BIDseq](https://github.com/y9c/pseudoU-BIDseq) 项目。

## 统计方法

统计分析方法统一了数据表示：对于每个候选 Ψ 位点，gap 校正后的 pileup 提供总深度（跨越位点的 reads 数）与缺失计数（包含缺失的 reads 数）。缺失计数占总深度的比例即为估计的 Ψ 比例，通过 beta-binomial 广义线性混合模型 (GLMM) 进行建模，以处理重复样本间的过离散（overdispersion）现象。

### Ψ 位点鉴定 ([`bin/Pseudo_callsites.R`](BID-transcript.script/bin/Pseudo_callsites.R))

**合并与预过滤**：样本对 (Mock/Treat) 按坐标（染色体、位置、链）合并。剔除组平均深度 < 10 的位点。

**组水平检验**：使用单侧 Fisher 精确检验 (Treat vs. Mock; 备择假设: "greater") 测试亚硫酸氢盐诱导的缺失富集。p 值经 Benjamini–Hochberg (BH) 校正；保留校正后 p < 0.01 的位点。

**重复标准**：位点需通过组水平检验，且在至少 `min_pass_reps` 个重复中满足以下条件：
- Treat 缺失计数 ≥ 5；
- Treat 缺失比例 ≥ 5%；
- Treat 缺失比例 ≥ 2 倍配对 Mock 重复。
*(对于重复数不等的样本，Treat 重复与 Mock 组平均值进行比较。)*

报告的比例与深度为通过筛选的重复项的平均值。

### PUS 依赖位点鉴定 ([`downstream/call_dependent_betabinomial.R`](BID-transcript.script/downstream/call_dependent_betabinomial.R))

**目的**：判断突变体中的修饰比例是否显著低于野生型 (WT)。

**流程**：
1. **候选池**：仅分析在 WT 中被判定为阳性的位点。
2. **深度阈值**：要求所有 WT 与突变体（Mock/Treat）重复的深度均 ≥ 10。
3. **统计模型**：使用 `glmmTMB` 拟合 beta-binomial GLMM。首先，通过位点水平随机截距估计全局离散度 (θ)，随后固定该参数。通过似然比 χ²(1) 检验 (`cbind(gap, depth − gap) ~ genotype`) 测试基因型效应。
4. **显著性**：基于 β 系数符号计算单侧 p 值，并进行 BH 校正；q < 0.05 且 Δ (WT_Ratio − Mut_Ratio) ≥ 效应量阈值（mRNA 与 tRNA 为 0.05，rRNA 为 0.10）的位点被视为依赖位点。
5. **最终判定**：仅保留在**两个**独立突变等位基因中均表现为依赖的位点。

### 温度响应位点鉴定 ([`downstream/call_temperature_response_betabinomial.R`](BID-transcript.script/downstream/call_temperature_response_betabinomial.R))

使用相同的 beta-binomial 框架，比较 WT 在 21 °C 与 27 °C 下的差异。候选池取两个温度下阳性位点的并集。执行双侧似然比检验：满足 BH q < 0.05 且 |Δ| ≥ 0.05 的位点判定为温度响应位点。

## 数据集与配置

[`configs/`](configs/) 包含以下五套数据集的 YAML 配置：

| 配置 | 流程 | 说明 |
|---|---|---|
| [`polyARNA.yaml`](configs/polyARNA.yaml) | polyA_mapping | poly(A)+ RNA：WT 与 7 个 PUS 基因的 14 个突变体株系（根、叶、花） |
| [`tRNA.yaml`](configs/tRNA.yaml) | totalRNA_mapping | tRNA 富集文库：基因型与组织同上 |
| [`Temperature.yaml`](configs/Temperature.yaml) | polyA_mapping | poly(A)+ RNA：WT 与 2 个 *trua1b* 株系（21 °C 与 27 °C） |
| [`Temperature_tRNA.yaml`](configs/Temperature_tRNA.yaml) | totalRNA_mapping | tRNA 富集文库：温度设计同上 |
| [`Chloroplast.yaml`](configs/Chloroplast.yaml) | totalRNA_mapping | 叶绿体 RNA |

### 配置参数说明
主要配置项：
- `paths`：辅助程序与比对索引的路径。
- `pipeline.ref`：参考序列 FASTA 文件（按比对顺序指定）。
- `pipeline.filter`：位点鉴定阈值（例如 `min_depth`, `max_padj`, `min_gap`, `min_pass_reps`）。
- `pipeline.umi_extract`：基于文库特异性 UMI/接头结构确定的 read 末端裁剪参数。
- `samples`：重复定义、FASTQ 路径及可选 `stages` (如 `[count]`) 以限定处理分支。

## 环境要求

### 分析流程
- **编排**：[Snakemake](https://snakemake.github.io)。
- **工具**：`trim_galore`, `fastp`, `seqkit`, `bowtie2`, `STAR`, `samtools`, `bedtools`, `csvtk`。
- **环境**：Python 3 (`polars`, `pysam`, `PyYAML`)。

### 下游分析
- R 包：`glmmTMB`, `optparse`, `data.table`, `jsonlite`。

## 数据可用性

原始测序数据将存入 accession **[TBD]**。处理后的位点特异性数据表位于 [`Figures/`](Figures/)。

## 引用
若您使用此代码或数据，请引用：
> [引用信息待定]
