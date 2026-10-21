# trua4_dependent_site_overview — WT vs trua4 叶绿体位点散点图

## 目的

只使用 leaf 组织，把叶绿体 rRNA / tRNA / mRNA 的全部 WT Ψ 位点画在
「X = trua4 突变体修饰比例，Y = WT 修饰比例」的散点图上，红点标记 `trua4` dependent 位点。
对应论文版 `trua4_dependent_site_overview`。

## 输入

WT leaf 位点：

- `polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv`
- `tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv`
- `polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv`

trua4 两个株系的修饰比例（`Treat_avg_ratio`）：

- `polyARNA/STAT/fisher/rRNA/trua4_{1,2}-leaf_sites.Fisher.tsv`
- `polyARNA/STAT/fisher/mRNA/trua4_{1,2}-leaf_sites.Fisher.tsv`
- `tRNA/STAT/fisher/ncRNA/trua4_{1,2}_Leaf_sites.Fisher.tsv`

dependent 标记（methylKit 口径，过滤 `is_dependent == TRUE`）：

- `polyARNA/STAT/dependent_methylKit/rRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
- `polyARNA/STAT/dependent_methylKit/mRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
- `tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`

叶绿体过滤：rRNA `Ch23S`、tRNA `^tRNA-Pt-`、mRNA `^ATCG`。

## 突变体比例合并规则

沿用论文版实现：

- 每个位点先取两个株系的 `Treat_avg_ratio`；
- 只把 `ratio > 0` 的株系计入，`Mut_Ratio` 为这些株系的均值；
- 两株系都为 0 时 `Mut_Ratio = 0`；
- 任一线系 flagged dependent，该位点即标为 `Dependent sites`
  （`dependency_support` 记录 `shared` / `trua4_1_only` / `trua4_2_only`）。

## 当前规模

| RNA | WT 位点 | dependent |
|---|---|---|
| rRNA | 4 | 2 |
| tRNA | 64 | 7 |
| mRNA | 18 | 0 |

## 输出

- `trua4_dependent_site_overview_dotplot.tsv`：每个位点的
  `WT_Ratio / Mut_Ratio / Delta_Ratio / site_status / dependency_support`
- `trua4_dependent_site_overview_summary.tsv`：每个 RNA 类型的位点数与 dependent 数
- `trua4_dependent_site_overview_dotplot_{rrna,trna,mrna}.pdf/png/svg`

散点图 0–100% 等比例坐标，虚线为 y = x；灰点 = WT sites（共享色板 `PAL_STATUS` 第一色），
红点 = Dependent sites（第二色）。

## 运行

```bash
bash data.sh
```
