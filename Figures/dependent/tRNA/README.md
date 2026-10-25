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
├── build_genotype_dependent_tables.sh    7 个基因型共用的数据构建逻辑（fix_pos 直接取上游 Sprinzl_pos）
├── dependent_trna_heatmap_common.R       热图公共实现（组织拆分 + 跨组织并集）
├── dependent_trna_distribution_common.R  Ψ 位点分布箱线图公共实现
├── sites_count/                          跨基因型：dependent 位点数堆叠柱状图
├── loss_effect/                          跨基因型：WT vs 突变体修饰比例散点
└── trua1a/ trua1b/ trua2/ trua3/ trua4/ trua5/ trua6/
    ├── data_prep.sh                      三行转发到 build_genotype_dependent_tables.sh
    ├── plot.R                            声明基因型 + 勾选本基因型要出的图
    └── <genotype>_*.{pdf,png,svg,tsv}    产物
```

代表 tRNA（每个 family 里 mean_CPM 最高的成员）来自
`tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`；本目录不再保存结构注释表。

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

- **`fix_pos` 直接取上游的 `Sprinzl_pos`**（canonical 0–76 位号）：上游
  `tRNA/STAT/dependent/ncRNA/<tissue>/<g>/<g>-<tissue>.dependent.Sprinzl.tsv` 已由
  `BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py` 基于 220 条完整注释生成，
  本目录不做任何坐标转换。上游按每个拷贝**自己的**区段注释计算，长度/区段划分与家族代表不同的拷贝也正确。
- dependent 表**没有**每株系修饰比 `<g>_1`/`<g>_2`，从
  `<g>_<line>-<tissue>.dependent-sites.tsv` 的 `Mut_Ratio` 取（该株系 3 个 Treat 重复的均值）。
- ncRNA 参考 `tRNA-unique_sn_snoRNA.fa` 混有 287 snoRNA + 82 snRNA，**必须按 `^tRNA-` 筛掉**。
- 热图行 = tRNA **家族**，代表成员取自 `family_top`；同家族所有拷贝的位点都并入该行。

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

- 行 = 组织 × tRNA 家族，按 Root / Leaf / Flowers 分三块（右侧色块标注组织）；
  家族 = 分室 + 氨基酸 + 反密码子，代表成员来自 `family_top`（mean_CPM 最高）
- 列 = canonical 0–76 位，按 p1 / Acceptor stem / D stem-loop / Anticodon stem-loop / Variable /
  T stem-loop / Acceptor stem / CCA 分块（底部有色带、位号、区段名）
- 值 = `Delta Ratio(WT - Mut)`，`#FFF5EE → #8B0000`；无位点的格子是浅灰
- 左侧窄色条 = 分室（Cy 橙 `#D55E00` / Mt 蓝 `#0072B2` / Pt 绿）；Cy/Mt/Pt 三分室全进图
- 图高随行数自适应 `max(8.8, 0.19*nrow+2.8)`；**同家族所有拷贝的位点都进图**，
  同一 (组织, 家族, canonical 位) 取 Delta 最大，不再丢弃非代表拷贝

**跨组织并集热图**（`trua1b` / `trua4` / `trua5`；`<g>_delta_ratio_union_dependent_trna_heatmap`）

- 行 = tRNA 家族，值 = 该位点跨三个组织**最大**的 Delta Ratio，便于看"最坏情况"

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

| 基因型 | all_tissues 行 | 并集位点 | 组织拆分 panel 行 | 并集家族数 | 图 |
|---|---|---|---|---|---|
| trua1a | 264 | 123 | 93 | 45 | 组织拆分 + cy_27/28 + cy_66/67/68 |
| trua1b | 79 | 60 | 51 | 38 | 组织拆分 + 并集 |
| trua2 | 87 | 56 | 74 | 42 | 组织拆分 + mt_27/28 |
| trua3 | 130 | 74 | 67 | 42 | 组织拆分 + cy_38/39/40 |
| trua4 | 60 | 47 | 47 | 34 | 组织拆分 + 并集 |
| trua5 | 100 | 66 | 56 | 40 | 组织拆分 + 并集 + cy_28 |
| trua6 | 232 | 99 | 89 | 42 | 组织拆分 + cy_38/39/40 + mt_38/39/40 |

旧版按“代表拷贝”过滤时，组织拆分 panel 行分别为 79 / 36 / 59 / 54 / 37 / 44 / 87；
新版把所有拷贝并回家族后增加到上表值（并集家族数同样从 37 / 27 / 29 / 33 / 26 / 31 / 40 增加）。

trua5 的 cy_28 面板：45 个 fix_pos==28 位点（Root 13 / leaf 15 / Flower 17），
效应量与 methylKit 版同量级（WT ~50% → 株系 ~13–28%）。

loss_effect 背景全量 = WT 阳性 tRNA 位点（Root 417 / leaf 505 / Flower 499），共 9947 行。

## 已知无害提示

- **`PostScript字体数据库里找不到'Arial'`**：`ComplexHeatmap::Heatmap()` 在"无设备"状态下测量文本宽度时，
  R 会临时打开默认 pdf 设备，其字体库里没有 Arial。测量值不受影响，
  且"构造一次、三种格式共用"保证 pdf/png/svg 布局完全一致。两个 common.R 里已定向屏蔽这一条，其余警告照常抛出。
- **`Using ragg device as default. Ignoring type and antialias`**：来自共享的 `save_plot()`，
  仓库所有图都有，PNG 实际走 ragg，Arial 经 fontconfig 正常解析。

## 验证记录（2026-10-25 重跑）

- 数据侧：7 个基因型的 `all_tissues` 行数与 `sites_count/dependent_site_count.tsv` 完全一致。
- 上游注释：`tRNA/STAT/dependent/ncRNA` 的 21 个最终 dependent 文件已生成 `.Sprinzl.tsv`，
  共 1034 行，其中 952 个 tRNA 行 100% 带 `Sprinzl_pos`，82 个 sn/snoRNA 行留空（被 `^tRNA-` 过滤）。
- 与旧布局（家族代表 + 银行家舍入）逐行对拍，差异行数：
  trua1a 0、trua1b 2、trua2 0、trua3 4、trua4 0、trua5 1、trua6 3；
  均为 `Asn-GTT-5`（57: 56→55）、`Leu-TAA-3`（40: 40→39）、`Ser-AGA-5`（41: 41→39）。
  其中 trua6 的 `Ser-AGA-5` 因此新进入 Ψ38/39/40 面板（+3 位点，n_values 540→549）。
- `fix_pos` 映射：`build_genotype_dependent_tables.sh` 对缺失/非法 `Sprinzl_pos` 直接 exit 1，
  7 个基因型全部零 miss 跑通。
- 图数：与 methylKit 版逐图对应（7 个基因型 17 张 + sites_count + loss_effect 21 张）。
