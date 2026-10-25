# merge_heatmap

## 目的

按组织分别整合 7 个基因型（trua1a/trua1b/trua2/trua3/trua4/trua5/trua6）的 dependent tRNA Ψ 位点，
以并集方式展示这些位点在 WT 和各突变体中的修饰比例。

移植自 `archive/20260902/paper/tRNA_dependent/merge_heatmap`，样式走 `Figures/shared`。

## 图形设计

- Y 轴：tRNA 家族（分室 + 氨基酸 + 反密码子），代表成员取自
  `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`（每个家族 mean_CPM 最高的 tRNA），按 Cy/Mt/Pt 分室分组；
  同家族所有拷贝的位点并入同一行，不再丢弃非代表拷贝。
- X 轴：校正后的 canonical 位号 `fix_pos`；每个位号下依次展示 WT、trua1a、trua1b、trua2、trua3、
  trua4、trua5、trua6 的修饰比例；位号标签格式 `tRNA_<fix_pos>`。
- 位点集合：该组织内所有突变体 `../<genotype>/<genotype>_all_tissues.dependent.tsv` 的并集；
  dependent 格标 `*`。
- 修饰比例：dependent 表只用于定义并集位点；WT 与突变体比例全部回填自
  `tRNA/STAT/adjustGap/ncRNA/` 的 `Treat-*.adjusted.tsv` 全位点 ratio 表（`data_prep.sh` 聚合）。
  因此即使某位点在某突变体中与 WT 差异不大，也显示该突变体的修饰比例。
- 颜色：修饰比例 0–1，缺失浅灰（`#FCFCFC`）。

## 与 archive 版的三处差异

| | archive | 本目录 |
|---|---|---|
| 位点来源 | Fisher 口径 dependent 表 | `Figures/dependent_betabinomial/tRNA/<g>/` 的 all_tissues 表（beta-binomial 口径，自带 fix_pos） |
| WT 比例 | `temp_totalRNA/STAT/fisher/.../WT-<tissue>_sites.Fisher.tsv` 的 Treat_avg_ratio | `Treat-WT-<tissue>.rep{1,2,3}.adjusted.tsv` 三重复均值——tRNA 流程的 Fisher.tsv 是过滤子集（Root 仅 623 位点）无法全位点回填；已抽查验证该均值与 Fisher.tsv 的 Treat_avg_ratio 逐位点一致 |
| 分室范围 | 只 Cy/Mt | Cy/Mt/Pt 全进图（2026-10-12 老板拍板） |

重复集口径与位点鉴定严格对齐：WT 的 rep4 是跨批次 count-only 重复，不进鉴定口径，
这里同样排除；突变体每株系 3 重复均值 → 两株系再均值，与 dependent 表
`<g>_1`/`<g>_2`/`Mut_Ratio` 的计算链一致（已抽查逐位点相等）。

## 输入

- `../trua*/trua*_all_tissues.dependent.tsv`（先跑各基因型目录的 `data_prep.sh`）
- `tRNA/STAT/adjustGap/ncRNA/Treat-{WT,<genotype>_<line>}-<tissue>.rep{1..3}.adjusted.tsv`
- `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`（家族代表，表达量口径）

## 输出

- `root/leaf/flower_merged_dependent_site_ratio_heatmap.{pdf,png,svg}`
- 对应的 `_long.tsv`（位点 × 样本长表）与 `_matrix.tsv`（比例矩阵）
- `merged_sample_ratios.tsv`（data_prep.sh 的聚合产物，plot.R 的唯一比例输入）

## 运行

```bash
bash data_prep.sh   # 聚合 adjusted.tsv → merged_sample_ratios.tsv（约 150 万行）
Rscript plot.R      # 画三张热图
```

规模参考（2026-10-25，beta-binomial 口径）：Root 517 并集位点（47 家族 × 11 位号）、leaf 819（63 × 13）、
Flower 696（58 × 12）；按家族合并后**不再丢弃任何拷贝**。
