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

**唯一上游**：`tRNA/STAT/dependent_methylKit/ncRNA/<tissue>/<genotype>/`
（由 `tRNA/STAT/dependent_methylKit/run.sh` 调用 `BID-transcript.script/downstream/call_dependent_methylKit.R` 生成）。

dependent 的判定口径（在 `call_dependent_methylKit.R` 里）：methylKit `F` 检验、单侧 p < 0.05、
效应量 `Delta_Ratio(WT-Mut) >= 0.05`、每株系 Mock 平均深度 >= 10，且**两个株系同时成立**（AND）。

### 与 archive 论文版的三处口径差异

| | archive（Fisher 版） | 本目录（methylKit 版） |
|---|---|---|
| dependent 表 | `dependent_site_overview/dependent/<g>-<tissue>.ncRNA.dependent.tsv` | `tRNA/STAT/dependent_methylKit/ncRNA/<tissue>/<g>/<g>-<tissue>.dependent.tsv` |
| canonical 位号 `fix_pos` | 表里自带 | **表里没有**，由 `data_prep.sh` 用 `All-tRNA_rep-structure.csv` 现算 |
| 每株系修饰比 `<g>_1`/`<g>_2` | 表里自带 | **表里没有**，从 `<g>_<line>-<tissue>.dependent-sites.tsv` 的 `Mut_Ratio` 取（该株系 3 个 Treat 重复的均值） |
| sn/snoRNA | 上游已按 `^tRNA` 过滤 | 参考 `tRNA-unique_sn_snoRNA.fa` 混有 287 snoRNA + 82 snRNA，**必须按 `^tRNA-` 筛掉** |

`fix_pos` = canonical 0–76 位号：按 `All-tRNA_rep-structure.csv` 的 17 个区段顺序与固定长度，
把 native pos 重标定为 `round((region_pos-1) * (fixed_len-1) / (region_len-1)) + 1`（半值取上）。
该结构表坐标 == 参考 fasta 的 1-based 坐标（含 CCA 尾），已在 `Figures/wildtype_landscape/representative_nuclear_trna_heatmap/README.md` 实测确认。

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

- 行 = 组织 × tRNA 代表，按 Root / Leaf / Flower 分三块（右侧色块标注组织）
- 列 = canonical 0–76 位，按 p1 / Acceptor stem / D stem-loop / Anticodon stem-loop / Variable /
  T stem-loop / Acceptor stem / CCA 分块（底部有色带、位号、区段名）
- 值 = `Delta Ratio(WT - Mut)`，`#FFF5EE → #8B0000`；无位点的格子是浅灰
- 左侧窄色条 = 分室（Cy 橙 `#D55E00` / Mt 蓝 `#0072B2`）
- 只画核与线粒体 tRNA。非代表拷贝（同一家族里与代表序列不一致的那些）不进图，
  运行时会打印 `Dropped N non-representative tRNA copies…` 列出被丢的 ID

**跨组织并集热图**（`trua1b` / `trua4` / `trua5`；`<g>_delta_ratio_union_dependent_trna_heatmap`）

- 行 = tRNA 代表，值 = 该位点跨三个组织**最大**的 Delta Ratio，便于看"最坏情况"
- 质体（Pt）tRNA 只进数据表、不进图：Pt 的环长分布与核/线粒体不同，投影到同一 canonical 布局有扭曲风险

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

| 基因型 | all_tissues 行 | 并集位点 | 热图 panel 行（组织 × 代表） | 图 |
|---|---|---|---|---|
| trua1a | 253 | 108 | 76 | 组织拆分 + cy_27/28 + cy_66/67/68 |
| trua1b | 24 | 19 | 10 | 组织拆分 + 并集 |
| trua2 | 49 | 23 | 25 | 组织拆分 + mt_27/28 |
| trua3 | 88 | 39 | 30 | 组织拆分 + cy_38/39/40 |
| trua4 | 24 | 13 | 2 | 组织拆分 + 并集（代表 tRNA 只有 1 个，图很矮） |
| trua5 | 61 | 31 | 19 | 组织拆分 + 并集 + cy_28 |
| trua6 | 209 | 82 | 68 | 组织拆分 + cy_38/39/40 + mt_38/39/40 |

## 已知无害提示

- **`PostScript字体数据库里找不到'Arial'`**：`ComplexHeatmap::Heatmap()` 在"无设备"状态下测量文本宽度时，
  R 会临时打开默认 pdf 设备，其字体库里没有 Arial。测量值不受影响（实测与 cairo 下一致），
  且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。两个 common.R 里已定向屏蔽这一条，其余警告照常抛出。
- **`Using ragg device as default. Ignoring type and antialias`**：来自共享的 `save_plot()`，
  仓库所有图都有，PNG 实际走 ragg，Arial 经 fontconfig 正常解析（已验证 PDF 内嵌 Arial-BoldMT、SVG 全 `font-family: "Arial"`、无 DejaVu）。
- **`Dropped N non-representative tRNA copies…`**：见上「组织拆分热图」，属预期行为。

## 验证记录

- 数据侧：7 个基因型的 `all_tissues` 行数与 `sites_count/dependent_site_count.tsv` 完全一致；
  并集行数 == 该表 distinct 位点数。
- 单点溯源：`trua1a` / `tRNA-Nu-Arg-CCG-1` / pos 68 逐值对上管线原表与两个株系的
  `dependent-sites.tsv`（`Mut_Ratio = (line1 + line2) / 2`）。
- `fix_pos` 映射：对全部 21 个「基因型×组织」文件的**每一个** tRNA 位点零 miss（7274 条映射）。
- 字体：新产 PDF 只含 Arial（Bold / Italic 子集），SVG 无 DejaVu；热图布局与 archive 版逐像素级一致。
