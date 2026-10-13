# dependent/tRNA — 按基因型的 tRNA dependent 位点图组（纯效应量口径）

## 目的

把每个 trua 突变体在 Root / Leaf / Flower 三个组织里**被该 TRUA 依赖**的 tRNA Ψ 位点摊开看：
位点落在 tRNA 二级结构的哪个位置、突变后掉多少、两个株系是否一致。

跨基因型的汇总（位点数 `sites_count`、修饰下降效应 `loss_effect`）在相邻目录；本目录是**按基因型展开**的 7 套图。

本目录与 `Figures/dependent_methylKit/tRNA/` 的图组**逐图同构**（脚本只差数据源指向），
数据源从 methylKit 口径换成 `tRNA/STAT/dependent`（纯效应量口径）。

## 目录结构

```
Figures/dependent/tRNA/
├── All-tRNA_rep-structure.csv            93 个 tRNA 代表的结构区段序列（46 Nu / 17 Mt / 30 Pt）
├── build_genotype_dependent_tables.sh    7 个基因型共用的数据构建逻辑
├── dependent_trna_heatmap_common.R       热图公共实现（组织拆分 + 跨组织并集）
├── dependent_trna_distribution_common.R  Ψ 位点分布箱线图公共实现
├── sites_count/                          跨基因型：dependent 位点数堆叠柱状图
├── loss_effect/                          跨基因型：WT vs 突变体修饰比例散点
└── trua1a/ trua1b/ trua2/ trua3/ trua4/ trua5/ trua6/
    ├── data_prep.sh                      三行转发到 build_genotype_dependent_tables.sh
    ├── plot.R                            声明基因型 + 勾选本基因型要出的图
    └── <genotype>_*.{pdf,png,svg,tsv}    产物
```

## 数据来源与口径

**唯一上游**：`tRNA/STAT/dependent/ncRNA/<tissue>/<genotype>/`
（由 `tRNA/STAT/dependent/run.sh` 调用 `BID-transcript.script/downstream/call_dependent.py` 生成，
解释器 `/home/HL/miniforge3/envs/DeepPSU/bin/python`）。

dependent 的判定口径（在 `call_dependent.py` 里，**无显著性检验**）：效应量
`Delta_Ratio(WT-Mut) >= 0.05`，且每株系 Mock 平均深度 >= 10，**两个株系同时成立**（AND）。

与 `dependent_methylKit` 版的关系（位点级已核验，见 `tRNA/STAT` 两侧）：methylKit 的
dependent 位点是本口径的**完全子集**（21/21 panel，无一例外）——甲基化检验只做"收紧"、
不引入新位点。本口径 tRNA-only 总量 952 个位点 vs methylKit 版 708。

### 与 methylKit 版图组的三处差异（其余逐字节同构）

| | Figures/dependent_methylKit/tRNA | 本目录 |
|---|---|---|
| dependent 表 | `tRNA/STAT/dependent_methylKit/ncRNA/…`（methylKit F 检验 + 效应量） | `tRNA/STAT/dependent/ncRNA/…`（纯效应量 delta>=0.05） |
| dependent 表列 | 多两列 p.value | 无 p.value 列（脚本按列名取数，不受影响） |
| 位点规模 | tRNA-only 708 | tRNA-only **952**（每基因型都更多，见下表） |

`fix_pos`、每株系修饰比、sn/snoRNA 过滤三件事与 methylKit 版完全一致：

- dependent 表**没有** `fix_pos`（canonical 0–76 位号），由 `data_prep.sh` 用
  `All-tRNA_rep-structure.csv` 现算：按 17 个区段顺序与固定长度，把 native pos 重标定为
  `round((region_pos-1) * (fixed_len-1) / (region_len-1)) + 1`（半值取上）。
  该结构表坐标 == 参考 fasta 的 1-based 坐标（含 CCA 尾）。
- dependent 表**没有**每株系修饰比 `<g>_1`/`<g>_2`，从
  `<g>_<line>-<tissue>.dependent-sites.tsv` 的 `Mut_Ratio` 取（该株系 3 个 Treat 重复的均值）。
- ncRNA 参考 `tRNA-unique_sn_snoRNA.fa` 混有 287 snoRNA + 82 snRNA，**必须按 `^tRNA-` 筛掉**。

## 运行

每个基因型目录独立：

```bash
cd trua1a
bash data_prep.sh                       # 只产数据，不画图
Rscript plot.R                          # 只画图
```

`data_prep.sh` 产出两张表：

- `<g>_all_tissues.dependent.tsv` — 三组织纵向拼接，`^tRNA-` only
  `chr, pos, fix_pos, strand, context, WT_Ratio, Mut_Ratio, Delta_Ratio(WT-Mut), <g>_1, <g>_2, Tissue`
- `<g>_union.dependent.tsv` — 跨组织并集，每位点取 Delta 最大的那次观测
  `…, Tissues, n_tissues, peak_tissue`

`plot.R` 只读这两张表，不碰 `tRNA/STAT/` 原始数据。

## 图怎么读

**组织拆分热图**（7 个基因型都有，`<g>_delta_ratio_tissue_split_dependent_trna_heatmap`）

- 行 = 组织 × tRNA 代表，按 Root / Leaf / Flowers 分三块（右侧色块标注组织）
- 列 = canonical 0–76 位，按 p1 / Acceptor stem / D stem-loop / Anticodon stem-loop / Variable /
  T stem-loop / Acceptor stem / CCA 分块（底部有色带、位号、区段名）
- 值 = `Delta Ratio(WT - Mut)`，`#FFF5EE → #8B0000`；无位点的格子是浅灰
- 左侧窄色条 = 分室（Cy 橙 `#D55E00` / Mt 蓝 `#0072B2` / Pt 绿）；Cy/Mt/Pt 三分室全进图
  （Pt tRNA 环长与核/线粒体不同、投影 canonical 布局有扭曲，按既定口径接受）
- 图高随行数自适应 `max(8.8, 0.19*nrow+2.8)`；非代表拷贝不进图，
  运行时打印 `Dropped N non-representative tRNA copies…` 列出被丢的 ID

**跨组织并集热图**（`trua1b` / `trua4` / `trua5`；`<g>_delta_ratio_union_dependent_trna_heatmap`）

- 行 = tRNA 代表，值 = 该位点跨三个组织**最大**的 Delta Ratio，便于看"最坏情况"

**Ψ 位点分布箱线图**（`trua1a` / `trua2` / `trua3` / `trua5` / `trua6`）

- x = 组织，y = Ψ 修饰比例 (%)，三组：WT（灰）、`<g>_1`（粉）、`<g>_2`（红）
- x 轴固定为 Roots / Leaves / Flowers 三档（`drop = FALSE`），某组织没有 dependent 位点时就留空
- 覆盖的区室×位点组合：

  | 基因型 | 面板 | 生物学含义 |
  |---|---|---|
  | trua1a | `cy_27_28`、`cy_66_67_68` | TRUA1A 在细胞质 tRNA 上的 Ψ27/28 与 Ψ66/67/68 |
  | trua2 | `mt_27_28` | TRUA2 的线粒体 Ψ27/28（组织受限） |
  | trua3 | `cy_38_39_40` | TRUA3 的根部特异细胞质 Ψ38/39/40 |
  | trua5 | `cy_28` | TRUA5 的细胞质 Ψ28（anticodon stem 5' 端第一位） |
  | trua6 | `cy_38_39_40`、`mt_38_39_40` | TRUA6 的双定位广谱 Ψ38/39/40 |

## 当前规模（2026-10-13 生成）

`all_tissues` 行数 = 该基因型三组织的 tRNA dependent 位点合计，与 `sites_count/dependent_site_count.tsv` 逐值一致。

| 基因型 | all_tissues 行 | 并集位点 | 组织拆分 panel 行 | 图 |
|---|---|---|---|---|
| trua1a | 264 | 123 | 79 | 组织拆分 + cy_27/28 + cy_66/67/68 |
| trua1b | 79 | 60 | 36 | 组织拆分 + 并集 |
| trua2 | 87 | 56 | 59 | 组织拆分 + mt_27/28 |
| trua3 | 130 | 74 | 54 | 组织拆分 + cy_38/39/40 |
| trua4 | 60 | 47 | 37 | 组织拆分 + 并集 |
| trua5 | 100 | 66 | 44 | 组织拆分 + 并集 + cy_28 |
| trua6 | 232 | 99 | 87 | 组织拆分 + cy_38/39/40 + mt_38/39/40 |

trua5 的 cy_28 面板：45 个 fix_pos==28 位点（Root 13 / leaf 15 / Flower 17），
效应量与 methylKit 版同量级（WT ~50% → 株系 ~13–28%）。

loss_effect 背景全量 = WT 阳性 tRNA 位点（Root 417 / leaf 505 / Flower 499），共 9947 行。

## 已知无害提示

- **`PostScript字体数据库里找不到'Arial'`**：`ComplexHeatmap::Heatmap()` 在"无设备"状态下测量文本宽度时，
  R 会临时打开默认 pdf 设备，其字体库里没有 Arial。测量值不受影响，
  且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。两个 common.R 里已定向屏蔽这一条，其余警告照常抛出。
- **`Using ragg device as default. Ignoring type and antialias`**：来自共享的 `save_plot()`，
  仓库所有图都有，PNG 实际走 ragg，Arial 经 fontconfig 正常解析。
- **`Dropped N non-representative tRNA copies…`**：见上「组织拆分热图」，属预期行为。

## 验证记录（2026-10-13）

- 数据侧：7 个基因型的 `all_tissues` 行数与 `sites_count/dependent_site_count.tsv` 完全一致；
  各基因型 loss_effect 表的 dependent 计数（71/98/95 等）与 sites_count 逐值一致。
- 上游时效：`tRNA/STAT/dependent` 产物为 2026-09-27 生成；用 2026-10-08 刷新后的
  adjustGap/fisher 输入对 trua1a-Flower 试跑 `call_dependent.py`，四个输出文件与现存版**逐字节一致**，
  证明上游重写为等值重写、现存结果无需重跑。
- `fix_pos` 映射：`build_genotype_dependent_tables.sh` 对映射失败直接 exit 1，
  7 个基因型全部零 miss 跑通。
- 图数：7 个基因型共 17 张 + sites_count 1 张 + loss_effect 21 张 = 39 图 ×3 格式，
  与 methylKit 版逐图对应。
