# dependent_methylKit/tRNA — 按基因型的 tRNA dependent 位点图组

## 目的

把每个 trua 突变体在 Root / Leaf / Flower 三个组织里**被该 TRUA 依赖**的 tRNA Ψ 位点摊开看：
位点落在 tRNA 二级结构的哪个位置、突变后掉多少、两个株系是否一致。

跨基因型的汇总（位点数 `sites_count`、修饰下降效应 `loss_effect`）在相邻目录；本目录是**按基因型展开**的 7 套图。

移植自论文版 `archive/20260902/paper/tRNA_dependent/{trua1a..trua6}`，数据源从旧 Fisher 口径换成
`tRNA/STAT/dependent_methylKit`（methylKit 口径），绘图样式改走 `Figures/shared`。

## 目录结构

```
Figures/dependent_methylKit/tRNA/
├── build_genotype_dependent_tables.sh    7 个基因型共用的数据构建逻辑（fix_pos 直接取上游 Sprinzl_pos）
├── dependent_trna_heatmap_common.R       热图公共实现（组织拆分 + 跨组织并集）
├── dependent_trna_distribution_common.R  Ψ 位点分布箱线图公共实现
├── sites_count/                          跨基因型：dependent 位点数堆叠柱状图
├── loss_effect/                          跨基因型：WT vs 突变体修饰比例散点
├── motif/                                跨基因型：各 基因型×组织 dependent 位点 5-nt 上下文 seqlogo（0 位点组合跳过）
├── merge_heatmap/                        跨基因型：每组织 dependent 位点并集 × (WT+7 突变体) 修饰比例热图
└── trua1a/ trua1b/ trua2/ trua3/ trua4/ trua5/ trua6/
    ├── data_prep.sh                      三行转发到 build_genotype_dependent_tables.sh
    ├── plot.R                            声明基因型 + 勾选本基因型要出的图
    └── <genotype>_*.{pdf,png,svg,tsv}    产物
```

代表 tRNA（每个 family 里 mean_CPM 最高的成员）来自
`tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`；本目录不再保存结构注释表。

## 数据来源与口径

**唯一上游**：`tRNA/STAT/dependent_methylKit/ncRNA/<tissue>/<genotype>/`
（由 `tRNA/STAT/dependent_methylKit/run.sh` 调用 `BID-transcript.script/downstream/call_dependent_methylKit.R` 生成）。

dependent 的判定口径（在 `call_dependent_methylKit.R` 里）：methylKit `F` 检验、
BH q < 0.05（方向为正，meth.diff > 0）、效应量 `Delta_Ratio(WT-Mut) >= 0.05`、
每株系 Mock 平均深度 >= 10，且**两个株系同时成立**（AND）。

### 与 archive 论文版的三处口径差异

| | archive（Fisher 版） | 本目录（methylKit 版） |
|---|---|---|
| dependent 表 | `dependent_site_overview/dependent/<g>-<tissue>.ncRNA.dependent.tsv` | `tRNA/STAT/dependent_methylKit/ncRNA/<tissue>/<g>/<g>-<tissue>.dependent.Sprinzl.tsv` |
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
只有 28（41 个，Root 12 / leaf 14 / Flower 15），故不写成 trua1a 那样的合并面板。
效应量很大：WT ~48–52% → 两个株系 ~13–28%，三个组织一致。

## 当前规模

`all_tissues` 行数 = 该基因型三组织的 tRNA dependent 位点合计，与 `sites_count/dependent_site_count.tsv` 逐值一致。

| 基因型 | all_tissues 行 | 并集位点 | 热图 panel 行（组织 × 家族） | 并集家族数 | 图 |
|---|---|---|---|---|---|
| trua1a | 217 | 94 | 76 | — | 组织拆分 + cy_27/28 + cy_66/67/68 |
| trua1b | 6 | 6 | 6 | 6 | 组织拆分 + 并集 |
| trua2 | 26 | 15 | 26 | — | 组织拆分 + mt_27/28 |
| trua3 | 42 | 25 | 21 | — | 组织拆分 + cy_38/39/40 |
| trua4 | 11 | 7 | 11 | 7 | 组织拆分 + 并集 |
| trua5 | 36 | 21 | 15 | 10 | 组织拆分 + 并集 + cy_28 |
| trua6 | 130 | 66 | 41 | — | 组织拆分 + cy_38/39/40 + mt_38/39/40 |

旧版按“代表拷贝”过滤，这些 panel 行数分别是 66 / 5 / 26 / 13 / 11 / 14 / 41；新版把所有拷贝并回家族后，
trua1a、trua1b、trua3、trua5 分别多出 10 / 1 / 8 / 1 行，trua2、trua4、trua6 不变。

## 已知无害提示

- **`PostScript字体数据库里找不到'Arial'`**：`ComplexHeatmap::Heatmap()` 在"无设备"状态下测量文本宽度时，
  R 会临时打开默认 pdf 设备，其字体库里没有 Arial。测量值不受影响（实测与 cairo 下一致），
  且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。两个 common.R 里已定向屏蔽这一条，其余警告照常抛出。
- **`Using ragg device as default. Ignoring type and antialias`**：来自共享的 `save_plot()`，
  仓库所有图都有，PNG 实际走 ragg，Arial 经 fontconfig 正常解析（已验证 PDF 内嵌 Arial-BoldMT、SVG 全 `font-family: "Arial"`、无 DejaVu）。

## 验证记录

- 数据侧：7 个基因型的 `all_tissues` 行数与 `sites_count/dependent_site_count.tsv` 完全一致；
  并集行数 == 该表 distinct 位点数。
- `fix_pos` 来源：21 个「基因型×组织」注释文件共 471 行，其中 468 个 tRNA 行
  **100% 带 `Sprinzl_pos`**（3 个 sn/snoRNA 行留空且被 `^tRNA-` 过滤）。
- 与旧版（家族代表布局）逐行对拍：464 行一致；4 行被修正——
  `trua3` 的 `Nu-Leu-TAA-3` pos 40 由 40 → 39（leaf、Flower），
  `trua6` 的 `Nu-Ser-AGA-5` pos 41 由 41 → 39（leaf、Flower）。
  后者因此新进入 trua6 的 Ψ38/39/40 面板（+2 位点，n_values 342 → 348）。
- 单点溯源：`trua1a` / `tRNA-Nu-Arg-CCG-1` / pos 68 逐值对上管线原表与两个株系的
  `dependent-sites.tsv`（`Mut_Ratio = (line1 + line2) / 2`）。
- 字体：新产 PDF 只含 Arial（Bold / Italic 子集），SVG 无 DejaVu。
