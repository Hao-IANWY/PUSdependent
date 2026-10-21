# chloroplast_purified_rna_validation — 叶绿体纯化 RNA 位点数量与 overlap 验证

## 目的

使用叶绿体纯化 RNA 建库数据，统计叶绿体 rRNA / tRNA / mRNA 上检测到的 Ψ 位点数量，
并与 WT leaf 的 poly(A)+ / small RNA 数据集做 overlap：
rRNA、mRNA 对 poly(A)+ 数据集，tRNA 对 small RNA 数据集。

对应论文版 `chloroplast_purified_rna_validation`。

## 输入

- `Chloroplast/STAT/fisher/ncRNA/Chloroplast_sites.positive.tsv`
- `Chloroplast/STAT/fisher/mRNA/Chloroplast_sites.positive.tsv`
- `Figures/chloroplast_RNA/chloroplast_site_overview/chloroplast_site_overview_density_sites.tsv`
  （WT leaf 对照集，由同级目录 `chloroplast_site_overview/data.sh` 生成）

## 过滤规则与 overlap 定义

- rRNA：从 ncRNA 结果中保留 `chr == "Ch23S"`
- tRNA：从 ncRNA 结果中保留 `chr` 以 `tRNA-Pt-` 开头
- mRNA：从 mRNA 结果中保留 `chr` 以 `ATCG` 开头
- 位点 ID：`chr:pos:strand`；overlap 为精确匹配

## 当前结果

| RNA | 纯化 RNA | WT leaf | overlap | 纯化 RNA 覆盖 | WT 覆盖 |
|---|---|---|---|---|---|
| rRNA | 3 | 4 | 3 | 100.0% | 75.0% |
| tRNA | 59 | 64 | 55 | 93.2% | 85.9% |
| mRNA | 15 | 18 | 11 | 73.3% | 61.1% |

## 输出

- `chloroplast_purified_rna_validation_sites.tsv`：纯化 RNA 过滤后的位点表
- `chloroplast_purified_rna_validation_count_summary.tsv`：位点数汇总
- `chloroplast_purified_rna_validation_overlap_summary.tsv`：overlap 汇总
  （含 `purified_only_sites` / `leaf_only_sites` / 双向覆盖率）
- `chloroplast_purified_rna_validation_barplot.pdf/png/svg`
- `chloroplast_purified_rna_validation_overlap.pdf/png/svg`：rRNA/tRNA/mRNA 三个 2-set Venn

Venn 用手绘 `ggplot` 圆（`geom_polygon`）实现，空集合也能正常出图。
柱状图的纵轴最高刻度由 `pretty()` 保证罩住最高柱，并额外留出柱顶数字空间。

## 运行

```bash
# 先跑 chloroplast_site_overview/data.sh 生成 WT leaf 对照表
bash data.sh
```
