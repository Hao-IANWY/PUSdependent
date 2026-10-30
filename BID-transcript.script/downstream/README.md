# downstream/ — 独立下游工具

本目录的工具**不被任何 smk 调用**，均手动运行。

## call_dependent_betabinomial.R — PUS 依赖位点判定（beta-binomial 版，主用）

PUS 依赖位点判定的主用统计实现。CLI 接口、深度门槛、效应量口径、
两株系 AND 判定与输出表结构（`<株系>.dependent-sites.tsv` / `<前缀>.dependent.tsv` /
`<前缀>.dependent.summary.json`）沿用原 methylKit 版的约定，可直接用同一套
run.sh 组织方式批量调用，检验用 Jie Huang（Beacon）
`dependency_stats_standalone.R`（见仓库 `dependent_jiehuang/`）的方法：

1. 在全部过深度门的位点上拟合含 site 随机截距的 beta-binomial GLMM，估全局
   离散度 theta；
2. 固定 theta 逐位点拟合 `cbind(gap, depth-gap) ~ genotype` 与 `~ 1` 两个
   模型，取 genotype 项的似然比 χ²(1) 双侧 P（`--site-test lrt` 默认，
   "WT 有、突变体全 0"的完全分离位点不退化；wald 为旧口径仅复现用），再按
   full 模型 β 符号折单侧 loss P（β<0 取 P/2）；
3. 单侧 P 在该株系全部候选位点内做 BH 校正，`q < --alpha` 且
   `Delta_Ratio(WT-Mut) >= --min-delta` 判 dependent。效应量门槛按 RNA
   类别区分：mRNA 0.05，rRNA 0.10（由各 `STAT/dependent_betabinomial/run.sh`
   控制，2026-10-30 起 rRNA 收紧）。

与原作者脚本的两处有意差异（保持本流程判定框架）：深度门槛沿用本流程严格
口径（突变株每个 Mock/Treat 重复、WT 每个重复都须 depth≥10；原方法只要求
每基因型 ≥2 个文库达标）；delta 用重复比例均值差（原方法用 pooled 差）。
逐位点 p/β 已对拍与原版 bit 级一致（tRNA leaf trua1a 567 候选、400 位点
抽检，theta 10 位小数相同）。

```bash
Rscript downstream/call_dependent_betabinomial.R \
    --tsv WT-Leaf_sites.Treat_avg_ratio.tsv \
    --key-cols chr,pos,strand --ratio-col Treat_avg_ratio \
    --genotype1 trua1a_1 --genotype2 trua1a_2 \
    --genotype1-treat-files Treat-trua1a_1.rep1.adjusted.tsv,...（rep1~rep3） \
    --genotype1-mock-files  Mock-trua1a_1.rep1.adjusted.tsv,...（同上） \
    --genotype2-treat-files Treat-trua1a_2.rep1.adjusted.tsv,...（同上） \
    --genotype2-mock-files  Mock-trua1a_2.rep1.adjusted.tsv,...（同上） \
    --outdir STAT/dependent_betabinomial/ncRNA/Leaf/trua1a \
    --output-prefix trua1a-Leaf
```

- `--tsv` 需含逐重复 Treat ratio/depth/**gap** 列（gap 用于计数；默认自动识别
  `Treat-*_gap`；Fisher 表原生带 `_gap` 列，无需预处理）。
- 默认 `--min-depth 10`、`--alpha 0.05`、`--min-delta 0.05`、`--significance q`、
  `--site-test lrt`、`--cores 1`；line 级表带 `beta_coef/beta_se/odds_ratio`
  审计列，summary JSON 记录两系 theta 与拟合失败计数。

## call_dependent_methylKit.R — PUS 依赖位点判定（methylKit 版，保留的参考实现）

比较 WT 与某 writer 基因的两个独立突变株系，判定 PUS 依赖位点（修饰在
突变体中显著丢失）。统计检验复用 methylKit 差异甲基化框架——把每个位点的
depth/gap 计数读作覆盖度/甲基化计数。**流程结果已迁移到上方 beta-binomial
版，项目内各数据集的 `STAT/dependent_methylKit/` 结果目录已删除**，本脚本
仅作为方法对照保留；原理与判定标准的完整描述见仓库顶层
[README 的「位点鉴定统计方法」](../../README.zh-CN.md#位点鉴定统计方法)一节。

```bash
Rscript downstream/call_dependent_methylKit.R \
    --tsv WT-Leaf_sites.Treat_avg_ratio.tsv \
    --key-cols chr,pos,strand --ratio-col Treat_avg_ratio \
    --genotype1 trua1a_1 --genotype2 trua1a_2 \
    --genotype1-treat-files Treat-trua1a_1.rep1.adjusted.tsv,Treat-trua1a_1.rep2.adjusted.tsv,Treat-trua1a_1.rep3.adjusted.tsv \
    --genotype1-mock-files  Mock-trua1a_1.rep1.adjusted.tsv,Mock-trua1a_1.rep2.adjusted.tsv,Mock-trua1a_1.rep3.adjusted.tsv \
    --genotype2-treat-files Treat-trua1a_2.rep1.adjusted.tsv,...（同上，rep1~rep3） \
    --genotype2-mock-files  Mock-trua1a_2.rep1.adjusted.tsv,...（同上，rep1~rep3） \
    --outdir STAT/dependent_methylKit/ncRNA/Leaf/trua1a \
    --output-prefix trua1a-Leaf
```

- `--tsv` 是 WT 候选位点表，须含逐重复 Treat ratio/depth 列（默认自动识别
  `Treat-*_ratio` / `Treat-*_depth`；列名不是该模式时须 `--wt-ratio-cols` /
  `--wt-depth-cols` 显式指定），因此不能用只有汇总列的 positive.tsv。
  项目内各数据集的 `STAT/dependent_methylKit/run.sh` 从 WT Fisher 全表抽取
  Pass==1 的行作此输入，并批量组织全部株系×组织的调用。
- 默认 `--min-depth 10`、`--alpha 0.05`、`--min-delta 0.05`、
  `--methylkit-test F`；输出逐株系 `<株系>.dependent-sites.tsv`、合并的
  `<前缀>.dependent.tsv` 与 `<前缀>.dependent.summary.json`。

## call_temperature_response_betabinomial.R — 温度响应位点判定（beta-binomial 版，主用）

比较同一 WT 材料在两个温度（如 27 °C vs 21 °C）下的修饰差异，判定温度
响应位点；检验框架与依赖版相同（beta-binomial GLMM + 似然比检验），但不
预设方向、直接使用双侧 P 做 BH 校正，方向（升/降）单列记录。候选位点为
两个温度 Fisher 阳性位点的并集（仅单一温度阳性的位点正是"另一温度修饰
丢失"的强候选），由各 `STAT/Temperature_response_betabinomial/run.sh` 合并
生成并批量调用；深度门槛、|Δ| ≥ 0.05 的效应量口径与依赖版一致。

```bash
Rscript downstream/call_temperature_response_betabinomial.R \
    --tsv WT_27_vs_WT_21.candidates.tsv \
    --key-cols chr,pos,strand \
    --group1 WT_21 --group2 WT_27 \
    --group1-treat-files Treat-21-WT.rep1.adjusted.tsv,Treat-21-WT.rep2.adjusted.tsv,Treat-21-WT.rep3.adjusted.tsv \
    --group1-mock-files  Mock-21-WT.rep1.adjusted.tsv,...（同上，rep1~rep3） \
    --group2-treat-files Treat-27-WT.rep1.adjusted.tsv,...（同上，rep1~rep3） \
    --group2-mock-files  Mock-27-WT.rep1.adjusted.tsv,...（同上，rep1~rep3） \
    --outdir STAT/Temperature_response_betabinomial/mRNA
```

- 输出 `<group2>_vs_<group1>.all-sites.tsv`、`.response.tsv`（含
  direction 列）与 `.response.summary.json`（含升/降计数）。

## call_dependent.py — PUS 依赖位点调用（纯效应量过渡版）

只用 Δratio 阈值判定、不做显著性检验的早期版本；正式分析使用上方带统计
检验的 `call_dependent_betabinomial.R`。以 WT 候选位点 TSV 为基准，比较两个
突变株系的 Treat/Mock adjusted TSV：

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

## fit_calibration.py — 校准曲线换算

给 TSV 的 ratio 列按 context 换算出 fraction 列（默认曲线为同目录的
`Calibration_curve.tsv`，可用 `--curve` 覆盖）：

```bash
python downstream/fit_calibration.py --tsv sites.Fisher.tsv --output sites.calibrated.tsv
```

## convert_to_transcript_coord.py — 基因组坐标 → 转录组坐标

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

### 使用要点

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

## extract_siteinfo.py — 提取指定位点的 depth/gap

`adjusted.tsv` 是全转录组逐碱基表（单文件四千万行以上），人工核对某几个位点时
不能整表打开。本工具按 (chr, pos) 精确匹配，把每个位点在多个样品上的
strand/context/ratio/depth/gap 汇总成一张表，查不到的组合记 NA。

输入是 yaml：`sites` 支持详细写法和 `"chr:pos"` 简写，`files` 支持 glob 通配符，
相对路径按 yaml 所在目录、仓库根目录两处依次解析：

```yaml
sites:
  - chr: AT1G06190.5          # 详细写法
    pos: 48
  - AT1G06190.5:66            # 简写
files:
  - Temperature/STAT/adjustGap/mRNA/Mock-21-trua1b_1.rep1.adjusted.tsv
  - Temperature/STAT/adjustGap/mRNA/Mock-27-*.rep?.adjusted.tsv   # glob
format: wide                  # 可选：wide（默认）/ long
output: siteinfo.tsv          # 可选；缺省打印到 stdout
```

```bash
python downstream/extract_siteinfo.py sites.yaml [-o out.tsv]
```

输出两种版式，`format` 控制：

- **wide（默认）** — 一个位点一行：`chr pos strand context` 之后每个样品展开
  `<sample>.ratio` / `<sample>.depth` / `<sample>.gap` 三列（N 样品即 4+3N 列）；
- **long** — 一个 site × sample 组合一行：`chr pos sample strand context ratio depth gap`。

样品名取文件名去掉 `.adjusted.tsv`，输入里出现同名样品时（跨 ref_type 混合）
用上级目录名做前缀消歧（`mRNA.X` / `rRNA.X`）。stderr 会报告每个文件的命中数。
文件可混跨 ref_type（mRNA/rRNA/ncRNA 的 chr 集不同，不在该文件里的位点记 NA）。
