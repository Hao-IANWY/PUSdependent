# dependent_betabinomial/tRNA/site_specficity — 位点特异性图组（按位号分文件夹）

## 目的

回答"WT 里在某个 canonical 位号（fixed 后的 Sprinzl 位置）上有修饰的**代表 tRNA**有多少个，
这些 tRNA 在一组目标突变体（各两条株系）里修饰比例变成什么样"。每个**位号**一个文件夹，
文件夹内固定两张图（一张热图 + 一张修饰比例分布图）。当前为 Flower 组的六个位号：

- `pos27/` `pos28/` `pos29/` — 基因型 trua1a + trua5（anticodon stem 区，1A 管 27/29、5 管 28）
- `pos38/` `pos39/` `pos40/` — 基因型 trua3 + trua6（anticodon loop 区，3 管 38、6 管 39/40）

**分室只保留核（Nu）**（老板 2026-10-28 拍板：忽略 Mt/Pt；图面不标分室，
分布图标题不带组织/n 标注）。

## 目录结构

```
site_specficity/
├── build_site_specificity_table.sh   共享数据构建（tissue × genotypes × sprinzl_pos × compartments 参数化）
├── site_specificity_common.R         热图 + 分布图公共实现（多基因型：每基因型两列比例 + 两列 per-line dependent）
├── WT-Flower_sites.Sprinzl.tsv       共享中间产物：WT positive 宇宙的 Sprinzl 注释版（同组织各位号文件夹共用）
├── pos27/ … pos40/
│   ├── data_prep.sh                  绑定 TISSUE=Flower / GENOTYPES / POS=N / COMPARTMENTS=Nu
│   ├── plot.R                        声明参数 + 勾选两张图
│   ├── pos<N>_flower_site_table.tsv  入选位点宽表（一行一个 tRNA 位点）
│   └── pos<N>_flower_site_specficity_{heatmap,distribution}.{pdf,png,svg}
```

加位号：复制任一文件夹，全局替换位号数字（必要时改 `GENOTYPES`）。加组织：复制文件夹改
`TISSUE` 与文件名（WT sprinzl 中间产物会按组织生成在共享层）。加基因型：两边同步加进
`GENOTYPES` / `genotypes`，并在 `site_specificity_common.R` 的
`site_specificity_line_palette` 补一对株系填充色（浅→深；现有 trua1a 粉/红、trua5
浅蓝/深蓝、trua3 浅绿/深绿、trua6 浅紫/深紫，缺色直接报错）。

## 入选口径

1. **WT 有修饰** = 位点在 WT 的 Fisher positive 宇宙
   （`tRNA/STAT/dependent_betabinomial/ncRNA/<tissue>/<g>/WT-<tissue>_sites.Treat_avg_ratio.tsv`，
   各基因型目录下逐字节相同，构建时 cmp 校验；表内全部行 Pass=1，即 WT 三个 Treat
   重复都检出 CMC 效应）。WT 比例取其 `Treat_avg_ratio`。
2. **位号** = `Sprinzl_pos == 目标位`（native pos 因 D 环/可变区长度不一而各异，
   如映射到 27 的 native pos 有 26/27/28/29）。WT 原表不带 Sprinzl 注释（上游只对
   21 张 final dependent 表做过），`data_prep.sh` 调 `add_sprinzl_pos.py` 在共享层补
   注释版中间产物。
3. **代表 tRNA** = `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`（每个 family
   mean_CPM 最高的成员）；`^tRNA-` 过滤 sn/snoRNA。
4. **分室** = `COMPARTMENTS` 绑定参数（Nu/Mt/Pt 逗号列表，空 = 全部）；当前只取
   `Nu`。Mt/Pt 不进表、不进图（Flower pos-27 全分室 28 个：Nu 13 / Mt 8 / Pt 7；
   其中 Mt/Pt 的 15 个在 trua1a 上全部不变，需要对照时可临时放开重跑）。

每条株系比例取 `<g>_<line>-<tissue>.dependent-sites.tsv` 的 `Mut_Ratio`（该株系 Treat
重复的均值，未筛选），与基因型图组一致；**株系表缺该位点时留 NA，绝不当作 0% 修饰**。
`data_prep.sh` 同时校验：行表 `WT_Ratio` 与 WT 宇宙 `Treat_avg_ratio`、两边 `context`
逐位点一致，否则报错退出。

## 图怎么读

**热图** `pos<N>_flower_site_specficity_heatmap`（样式对齐 `../../mRNA/heatmap`）

- 行 = 该位号 WT 修饰的核代表 tRNA，按 WT 比例降序；行标签 = family
  （`tRNA-Nu-` 前缀与拷贝号去掉，如 `Arg-CCT`）
- 列 = WT + 各基因型两条株系；格值 = 修饰比例（<10% 一位小数，其余取整）
- 红框 = 该株系 beta-binomial **per-line** 判定 dependent（单边 TRUE 也画框，与下面
  数字表"双株系 TRUE"的计数可以差 1）；灰格 NA = 该株系无测量

**修饰比例分布** `pos<N>_flower_site_specficity_distribution`
（样式对齐 `../trua1a` 的 position_distribution 箱线 + 抖动点）

- x = 五组样本，y = Ψ 修饰比例 (%)，单面板，标题只写 "Ψ<N> WT-modified tRNAs"
- WT 灰 / 各基因型浅→深一对色（trua1a 粉红、trua5 蓝、trua3 绿、trua6 紫）

## 当前数字（Flower，Nu-only）

| 位号 | tRNA 数 | trua1a | trua5 | trua3 | trua6 | 中位 WT → 塌方 |
|---|---|---|---|---|---|---|
| Ψ27 | 13 | **11** | 0 | — | — | 46.9% → 5.6/4.9% (1a) |
| Ψ28 | 4 | 1 | **4** | — | — | 40.6% → 12.8/18.4% (5) |
| Ψ29 | 1（Val-TAC） | **1** | 0 | — | — | 43.4% → 4.6/4.1% (1a) |
| Ψ38 | 5 | — | — | **5** | 0 | 26.9% → 7.3/4.7% (3) |
| Ψ39 | 16 | — | — | 1 | **12** | 61.8% → 13.3/14.9% (6) |
| Ψ40 | 1（Ser-GCT） | — | — | 0 | **1** | 73.9% → 15.0/18.8% (6) |

（加粗 = 双株系 simultaneously dependent 的位点数；"—" = 该位号图组未含此基因型）

两组互为镜像的特异性：Ψ27/29 归 TRUA1A、Ψ28 归 TRUA5（anticodon stem）；Ψ38 归
TRUA3、Ψ39/40 归 TRUA6（anticodon loop）。值得注意的单点：
- tRNA-Nu-Cys-GCA-10 的 Ψ27 在 trua1a 不掉，但 Ψ28 在 trua1a 和 trua5 双 dependent
  ——同一条 tRNA 上两个位号各归其主还有交叉；
- Ψ27 的 trua1a 例外 Arg-TCG-2（0.23 → 0.06，比例掉但 per-line q 未过门槛）；
- Ψ39 的 Trp-CCA-1 只在 trua6_1 单边过门槛（_2 FALSE），红框 27 vs 双 TRUE 26 的差。

## 运行

```bash
cd pos27          # 或其余位号目录
bash data_prep.sh # 只产数据，不画图
Rscript plot.R    # 只画图
```

## 验证记录

- 2026-10-28：pos27 的 4 条株系列（trua1a_1/2、trua5_1/2 比例 + per-line dependent）
  与上游 4 份 `dependent-sites.tsv` 逐位点对拍全 OK；11 个双株系 TRUE 的 trua1a 位点
  全部以 `Sprinzl_pos=27` 出现在 `trua1a-Flower.dependent.Sprinzl.tsv`（final 表另有
  3 行 Ψ65/67/68 属同批 tRNA，与本图口径无关）。
- 热图红框 SVG 计数核对：pos27 = 22、pos28 = 10、pos29 = 2、pos38 = 10、pos39 = 27
  （含 Trp-CCA-1 的 trua6_1 单边框）、pos40 = 2；均与 per-line TRUE 计数吻合。
- 与 `tRNA/STAT/fisher/ncRNA/WT_Flower_sites.positive.Sprinzl.tsv`（上游现成注释版，
  仅 ratio/depth 列）独立对数：pos-27 全分室 57 行、代表 tRNA 28（Nu 13），一致。
- 株系列当前无 NA；缺位点时热图灰格 + "NA"、分布图丢点并记数。
- 踩坑记录：awk END 拼 (值,dependent) 交替而表头是"两值在前两 dep 在后"，列错位导致
  readr 把比例列读成 logical/数值混杂——已改为与表头一致的双循环拼行。
