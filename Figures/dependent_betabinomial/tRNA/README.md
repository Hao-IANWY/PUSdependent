# dependent_betabinomial/tRNA — 按基因型的 tRNA dependent 位点图组（beta-binomial 口径）

## 目的

把每个 trua 突变体在 Root / Leaf / Flower 三个组织里**被该 TRUA 依赖**的 tRNA Ψ 位点摊开看：
位点落在 tRNA 二级结构的哪个位置、突变后掉多少、两个株系是否一致。

跨基因型的汇总（位点数 `sites_count`、修饰下降效应 `loss_effect`）在相邻目录；本目录是**按基因型展开**的 7 套图。

移植自论文版 `archive/20260902/paper/tRNA_dependent/{trua1a..trua6}`，数据源从旧 Fisher 口径换成
`tRNA/STAT/dependent_betabinomial`（beta-binomial 口径），绘图样式改走 `Figures/shared`。

## 目录结构

```
Figures/dependent_betabinomial/tRNA/
├── build_genotype_dependent_tables.sh    7 个基因型共用的数据构建逻辑（fix_pos 直接取上游 Sprinzl_pos）
├── dependent_trna_heatmap_common.R       热图公共实现（组织拆分 + 跨组织并集）
├── dependent_trna_distribution_common.R  Ψ 位点分布箱线图公共实现
├── dependent_trna_motif_common.R         位点 motif seqlogo 公共实现（与分布图同套选点口径）
├── sites_count/                          跨基因型：dependent 位点数堆叠柱状图
├── loss_effect/                          跨基因型：WT vs 突变体修饰比例散点
├── motif/                                跨基因型：各 基因型×组织 dependent 位点 5-nt 上下文 seqlogo（0 位点组合跳过）
├── merge_heatmap/                        跨基因型：每组织 dependent 位点并集 × (WT+7 突变体) 修饰比例热图
├── site_specficity/                      按位号分文件夹：WT 在该位号有修饰的核代表 tRNA × 对应突变体组（见其 README）
└── trua1a/ trua1b/ trua2/ trua3/ trua4/ trua5/ trua6/
    ├── data_prep.sh                      三行转发到 build_genotype_dependent_tables.sh
    ├── plot.R                            声明基因型 + 勾选本基因型要出的图
    └── <genotype>_*.{pdf,png,svg,tsv}    产物
```

代表 tRNA（每个 family 里 mean_CPM 最高的成员）来自
`tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`；本目录不再保存结构注释表。

## 数据来源与口径

**唯一上游**：`tRNA/STAT/dependent_betabinomial/ncRNA/<tissue>/<genotype>/`
（由 `tRNA/STAT/dependent_betabinomial/run.sh` 调用 `BID-transcript.script/downstream/call_dependent_betabinomial.R` 生成）。

dependent 的判定口径（在 `call_dependent_betabinomial.R` 里，Jie Huang beta-binomial
方法）：全局 theta（site 随机截距 beta-binomial GLMM）固定后逐位点 Wald 单侧 P、
BH q < 0.05、效应量 `Delta_Ratio(WT-Mut) >= 0.05`、每株系 Mock 平均深度 >= 10，
且**两个株系同时成立**（AND）。与 methylKit 版的对照图见 `sites_count/overlap/`。

### 与 archive 论文版的三处口径差异

| | archive（Fisher 版） | 本目录（beta-binomial 版） |
|---|---|---|
| dependent 表 | `dependent_site_overview/dependent/<g>-<tissue>.ncRNA.dependent.tsv` | `tRNA/STAT/dependent_betabinomial/ncRNA/<tissue>/<g>/<g>-<tissue>.dependent.Sprinzl.tsv` |
| canonical 位号 `fix_pos` | 表里自带 | **直接取表里的 `Sprinzl_pos`**（由 `BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py` 基于 220 条完整注释生成），本目录不做任何坐标转换 |
| 代表 tRNA | 未定义 | `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`（每个 family mean_CPM 最高的成员） |
| 每株系修饰比 `<g>_1`/`<g>_2` | 表里自带 | **表里没有**，从 `<g>_<line>-<tissue>.dependent-sites.tsv` 的 `Mut_Ratio` 取（该株系 3 个 Treat 重复的均值） |
| sn/snoRNA | 上游已按 `^tRNA` 过滤 | 参考 `tRNA-unique_sn_snoRNA.fa` 混有 287 snoRNA + 82 snRNA，**必须按 `^tRNA-` 筛掉** |

`fix_pos` 就是上游的 `Sprinzl_pos`（canonical 0–76，Sprinzl 编号）；上游按每个拷贝**自己的**区段注释计算，
因此长度/区段划分与家族代表不同的拷贝也正确（旧版会把这些拷贝的位点算偏或剔除）。

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

- 行 = 组织 × tRNA 家族，按 Root / Leaf / Flower 分三块（右侧色块标注组织）；
  家族 = 分室 + 氨基酸 + 反密码子，代表成员来自 `family_top`（mean_CPM 最高）
- 列 = canonical 0–76 位，按 p1 / Acceptor stem / D stem-loop / Anticodon stem-loop / Variable /
  T stem-loop / Acceptor stem / CCA 分块（底部有色带、位号、区段名）
- 值 = `Delta Ratio(WT - Mut)`，`#FFF5EE → #8B0000`；无位点的格子是浅灰
- 左侧窄色条 = 分室（Cy 橙 `#D55E00` / Mt 蓝 `#0072B2` / Pt 绿 `#009E73`）
- **同家族所有拷贝的位点都进图**：位点按 `chr` 归到家族，同一 (组织, 家族, canonical 位)
  收到多个拷贝时取 Delta 最大；不再丢弃非代表拷贝

**跨组织并集热图**（`trua1b` / `trua4` / `trua5`；`<g>_delta_ratio_union_dependent_trna_heatmap`）

- 行 = tRNA 家族，值 = 该位点跨三个组织**最大**的 Delta Ratio，便于看"最坏情况"
- Cy / Mt / Pt 三分室都进图（2026-10-12 老板拍板）

**Ψ 位点分布箱线图**（`trua1a` / `trua2` / `trua3` / `trua5` / `trua6`）

- x = 组织，y = Ψ 修饰比例 (%)，三组：WT（灰）、`<g>_1`（粉）、`<g>_2`（红）
- 每个点是一个位点，箱线图为该组分布；点是同一批位点在三个样本里的修饰比例
- x 轴固定为 Roots / Leaves / Flowers 三档（`drop = FALSE`），某组织没有 dependent 位点时就留空
- 覆盖的区室×位点组合：

  | 基因型 | 面板 | 生物学含义 |
  |---|---|---|
  | trua1a | `cy_27_28`、`cy_66_67_68` | TRUA1A 在细胞质 tRNA 上的 Ψ27/28 与 Ψ66/67/68 |
  | trua2 | `mt_27_28` | TRUA2 的线粒体 Ψ27/28（组织受限） |
  | trua3 | `cy_38_39_40` | TRUA3 的根部特异细胞质 Ψ38/39/40 |
  | trua5 | `cy_28` | TRUA5 的细胞质 Ψ28（anticodon stem 5' 端第一位） |
  | trua6 | `cy_38_39_40`、`mt_38_39_40` | TRUA6 的双定位广谱 Ψ38/39/40 |

**trua5 的 Ψ28 为什么是单点而不是 Ψ27/28**：trua5 在 canonical 27 上没有任何 dependent 位点，
只有 28（41 个，Root 10 / leaf 15 / Flower 16），故不写成 trua1a 那样的合并面板。
效应量很大：WT ~48–52% → 两个株系 ~13–28%，三个组织一致。

**位点 motif seqlogo**（`trua1a` / `trua2` / `trua3` / `trua5` / `trua6`，`<g>_<suffix>_dependent_motif`）

- 输入 `<g>_union.dependent.tsv`：跨组织去重，同一位点不因多组织重复而抬高 bits；
  上下文取 5-nt `context`，x 轴 -2..+2、中心 0 位 = 修饰位点（Ψ 来源 U，图中恒为 T）
- 选点口径与分布箱线图同构（compartment × Sprinzl 位号），每个组合一张 logo，样式与
  `motif/`（基因型×组织全景版）一致；位点数落盘 `<g>_dependent_motif_site_counts.tsv`
- 覆盖的组合（2026-10-27 老板指定；trua3 的 motif 只写 38/39，数据里本就没有 Cy_40
  位点，与其分布图的 38/39/40 选择等价）：

  | 基因型 | 面板 | 独立位点数 |
  |---|---|---|
  | trua1a | `cy_27_28`、`cy_66_67_68` | 41、19 |
  | trua2 | `mt_27_28` | 9 |
  | trua3 | `cy_38_39` | 22 |
  | trua5 | `cy_28` | 15 |
  | trua6 | `cy_38_39_40`、`mt_38_39_40` | 58、7（数据中实际落在 39/40） |

## 当前规模

`all_tissues` 行数 = 该基因型三组织的 tRNA dependent 位点合计，与 `sites_count/dependent_site_count.tsv` 逐值一致。

| 基因型 | all_tissues 行 | 并集位点 | 热图 panel 行（组织 × 家族） | 并集家族数 | 图 |
|---|---|---|---|---|---|
| trua1a | 179 | 72 | 54 | — | 组织拆分 + cy_27/28 + cy_66/67/68 + motif×2 |
| trua1b | 13 | 9 | 6 | 5 | 组织拆分 + 并集 |
| trua2 | 41 | 17 | 41 | — | 组织拆分 + mt_27/28 + motif×1 |
| trua3 | 62 | 25 | 28 | — | 组织拆分 + cy_38/39/40 + motif×1 |
| trua4 | 22 | 11 | 22 | 11 | 组织拆分 + 并集 |
| trua5 | 44 | 18 | 14 | 5 | 组织拆分 + 并集 + cy_28 + motif×1 |
| trua6 | 183 | 64 | 59 | — | 组织拆分 + cy_38/39/40 + mt_38/39/40 + motif×2 |

与 methylKit 版（`Figures/dependent_methylKit/tRNA`）的逐组织计数对照见
`sites_count/overlap/dependent_site_overlap_summary.tsv`：tRNA-only 合计
beta-binomial 544 / methylKit 468 / 交集 397。最大差异在 Root（154 vs 44）：
methylKit 在 Root 的 line2 过离散 φ 爆炸把功效吃光（同表大量位点 q 齐卡
0.055），beta-binomial 的全局 theta 保留了功效（详见上游
`tRNA/STAT/dependent_betabinomial` 的说明）。

## 已知无害提示

- **`PostScript字体数据库里找不到'Arial'`**：`ComplexHeatmap::Heatmap()` 在"无设备"状态下测量文本宽度时，
  R 会临时打开默认 pdf 设备，其字体库里没有 Arial。测量值不受影响（实测与 cairo 下一致），
  且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。两个 common.R 里已定向屏蔽这一条，其余警告照常抛出。
- **`Using ragg device as default. Ignoring type and antialias`**：来自共享的 `save_plot()`，
  仓库所有图都有，PNG 实际走 ragg，Arial 经 fontconfig 正常解析（已验证 PDF 内嵌 Arial-BoldMT、SVG 全 `font-family: "Arial"`、无 DejaVu）。

## 验证记录

- motif seqlogo：`<g>_dependent_motif_site_counts.tsv` 的 n 与 union 表按
  compartment × fix_pos 手工筛选一致（41/19/9/22/15/58/7），各组合中心碱基 100% 为
  T（Ψ 来源 U）；运行会连带重刷同目录热图与分布图（同输入确定性重生成，内容不变）。
- 数据侧：7 个基因型的 `all_tissues` 行数与 `sites_count/dependent_site_count.tsv`
  完全一致（179/13/41/62/22/44/183，合计 544 = STAT summary 547 − 3 个 sn/sno
  位点，`^tRNA-` 过滤口径）；并集行数 == 该表 distinct 位点数。
- `fix_pos` 来源：Sprinzl 注入由 `add_sprinzl_pos.py` 对
  `tRNA/STAT/dependent_betabinomial` 的 21 张 final 表逐张执行，映射规则与
  methylKit 版共用同一套按拷贝注释。
- overlap Venn 的两方法集合 = `STAT/dependent_betabinomial` vs
  `STAT/dependent`（纯效应量版，与 methylKit 图组的对比对象同构）；
  beta-binomial 位点为其完全子集（tRNA-only 544 ⊂ 952）。
- 字体与布局沿用 `Figures/shared`；新产 PDF 只含 Arial 子集。
