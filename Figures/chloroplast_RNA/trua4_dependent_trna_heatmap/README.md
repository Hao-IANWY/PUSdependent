# trua4_dependent_trna_heatmap — Leaf 叶绿体 dependent tRNA 结构热图

## 目的

只使用 Leaf 组织中 `trua4` dependent 的叶绿体 tRNA 位点，展示这些位点在 canonical
tRNA 结构（0–76 位）上的分布与修饰下降幅度。对应论文版
`archive/20260902/paper/chloroplast_RNA/trua4_dependent_trna_heatmap`。

## 输入

数据源为 beta-binomial 口径（`STAT/dependent_betabinomial`）：

- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4-leaf.dependent.Sprinzl.tsv`
  （官方 dependent 表，两株系 AND；带 `Sprinzl_pos`）
- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
  （两株系 `Mut_Ratio`，只作列记录）

## 过滤与映射规则

- 只保留 `chr` 以 `tRNA-Pt-` 开头的叶绿体 tRNA（排除同一 ncRNA 参考里的核/线粒体
  tRNA 与 sn/snoRNA）
- canonical 位号直接取表里的 `Sprinzl_pos`，不再依赖旧的 `All-tRNA_rep-structure.csv`
- 列布局 = canonical 0–76 位：p1 / Acceptor stem / D stem-loop / Anticodon stem-loop /
  Variable / T stem-loop / Acceptor stem / CCA，共 17 个区段、77 列，与
  `Figures/dependent_betabinomial/tRNA` 的热图同构
- 行 = 每个 dependent 叶绿体 tRNA；按 canonical 位号升序、再按 tRNA 名称排序，
  使 Ψ38 / Ψ39 / Ψ40 三块一目了然
- 值 = `Delta_Ratio(WT-Mut)`；同一 (tRNA, canonical 位) 有多个观测时取最大值

## 当前规模

9 个 dependent 位点，全部为两株系共享（`shared`）：

| canonical 位号 | tRNA |
|---|---|
| 38 | Pt-Pro-TGG-1 |
| 39 | Pt-Cys-GCA-1、Pt-His-GTG-1、Pt-Leu-CAA-1、Pt-Leu-TAA-1、Pt-Ser-GCT-1、Pt-Thr-TGT-1、Pt-Tyr-GTA-1 |
| 40 | Pt-Asn-GTT-1 |

即 TRUA4 在叶绿体 tRNA 上只作用于反密码子环 3' 端 / 反密码子茎的 Ψ38–40。

与论文版热图的差异：论文版是 7 个位点（6 个共享 + Pt-Val-TAC Ψ55 的
`trua4_1_only`）；beta-binomial 口径下 Val-TAC 两个株系均不 dependent
（q = 0.63 / 0.87），而 Cys-GCA、Leu-TAA、Ser-GCT 三个位点新进入 dependent，
故当前为 9 个位点、全部落在 Ψ38–40，核心结论不变。

## 输出

- `trua4_pt_trna_dependent_union.tsv`：叶绿体 dependent 位点表
  （`chr, pos, Sprinzl_pos, strand, context, WT_Ratio, Mut_Ratio,
  Delta_Ratio(WT-Mut), trua4_1, trua4_2, dependency_support`）
- `trua4_pt_trna_dependent_mapped_sites.tsv`：叠加 `display_label`、`canonical_pos`、
  `matrix_col_id` 后的映射表
- `trua4_pt_trna_dependent_heatmap.pdf/png/svg`

图注：行为每个 dependent 叶绿体 tRNA，列为 canonical 0–76 位；颜色 = `Delta Ratio`
（`#FFF5EE → #8B0000`），色标上限由数据最大值经 `y_axis_cover()` 取整，保证覆盖最大值；
无位点的格子为浅灰。

## 运行

```bash
bash data.sh
```
