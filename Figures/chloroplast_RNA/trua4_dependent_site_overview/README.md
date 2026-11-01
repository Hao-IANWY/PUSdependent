# trua4_dependent_site_overview — Leaf 叶绿体位点 WT vs trua4 散点

## 目的

只使用 Leaf 组织中的 WT 与 `trua4` 数据，展示全部叶绿体 rRNA / tRNA / mRNA 位点在
WT 和 `trua4` 中的修饰比例，并用红点标记 `trua4` dependent 位点。对应论文版
`archive/20260902/paper/chloroplast_RNA/trua4_dependent_site_overview`。

## 输入

数据源为 beta-binomial 口径（`STAT/dependent_betabinomial`），与
`Figures/dependent_betabinomial` 一致：

- `polyARNA/STAT/dependent_betabinomial/rRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv`
- `polyARNA/STAT/dependent_betabinomial/rRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv`
- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
- `polyARNA/STAT/dependent_betabinomial/mRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv`
- `polyARNA/STAT/dependent_betabinomial/mRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`

WT 表就是 WT 阳性位点表（`WT-*_sites.positive.tsv` 的 `Treat_avg_ratio` 版本，
位点集与 `ratio` 逐值相同）。

## 过滤规则

- rRNA：`chr == "Ch23S"`
- tRNA：`chr` 以 `tRNA-Pt-` 开头（排除同一 ncRNA 参考里的核/线粒体 tRNA 与 sn/snoRNA）
- mRNA：`chr` 以 `ATCG` 开头（叶绿体编码基因）

## dependent 口径

与 `Figures/dependent_betabinomial` 一致：**两个株系 `is_dependent` 同时为 TRUE**
才标红（`site_status = Dependent sites`）。只在一株系成立的位点在
`dependency_support` 里记为 `trua4_1_only` / `trua4_2_only`，但不上色。
`Mut_Ratio` = 两株系 `Mut_Ratio` 的均值。

## 当前规模

| RNA | WT 位点 | dependent |
|---|---|---|
| rRNA | 4 | 2 |
| tRNA | 64 | 9 |
| mRNA | 18 | 0 |

tRNA 另有 1 个位点只在 `trua4_2` 成立（`trua4_2_only`），按 AND 口径不计入 dependent。

## 输出

- `trua4_dependent_site_overview_dotplot.tsv`：三个 RNA 类型合并的位点表
  （`rna_type, chr, pos, strand, site_id, WT_Ratio, Mut_Ratio, Delta_Ratio,
  site_status, dependency_support, trua4_1, trua4_2, trua4_1_dep, trua4_2_dep`）
- `trua4_dependent_site_overview_summary.tsv`：每个 RNA 类型的位点数与 dependent 数
- `trua4_dependent_site_overview_dotplot_{rrna,trna,mrna}.pdf/png/svg`

图注：横轴 = `trua4` 突变体平均修饰比例，纵轴 = WT 平均修饰比例；虚线为 y = x；
灰点 = 全部 WT 位点，红点 = dependent 位点；坐标轴 0–100%。

## 运行

```bash
bash data.sh
```
