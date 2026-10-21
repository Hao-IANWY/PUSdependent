# chloroplast_site_overview — WT leaf 叶绿体 Ψ 位点总览

## 目的

只使用 WT leaf 组织，汇总叶绿体 rRNA / tRNA / mRNA 上全部 Ψ 位点的数量与修饰比例分布，
每个 RNA 类型一张密度图。对应论文版 `chloroplast_site_overview`。

## 输入

- `polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv`
- `tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv`
- `polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv`

## 过滤规则

- rRNA：`chr == "Ch23S"`
- tRNA：`chr` 以 `tRNA-Pt-` 开头
- mRNA：`chr` 以 `ATCG` 开头
- 保留 `ratio` 为有限值且落在 0–1 的位点，展示时换算为百分比

## 当前规模

| RNA | 位点数 |
|---|---|
| rRNA | 4 |
| tRNA | 64 |
| mRNA | 18 |

## 输出

- `chloroplast_site_overview_density_sites.tsv`：三个 RNA 类型合并的位点表
  （`rna_type, chr, pos, strand, ratio, depth, ratio_pct, site_id`）
- `chloroplast_site_overview_count_summary.tsv`：每个 RNA 类型的位点数
- `chloroplast_site_overview_density_{rrna,trna,mrna}.pdf/png/svg`

图注：横轴 = 修饰比例（0–100%），纵轴 = 密度，曲线内标注 `n = 位点数`。
纵轴刻度由 `y_axis_cover()` 生成，最高刻度保证罩住密度峰值。

## 运行

```bash
bash data.sh
```
