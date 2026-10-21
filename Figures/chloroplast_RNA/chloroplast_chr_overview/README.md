# chloroplast_chr_overview — 叶绿体染色体（ChrC）环形 Ψ 位点总览

## 目的

基于 `Araport11_fixed.gtf` 的 `ChrC` 注释绘制拟南芥叶绿体染色体环形图：
外圈标注叶绿体 rRNA / tRNA / mRNA 基因位置，内圈用径向 bar 展示 WT leaf 中 Ψ 位点的
修饰比例。对应论文版 `chloroplast_chr_overview`。

## 输入

- `Araport11_fixed.gtf`
- `polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv`
- `tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv`
- `polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv`

## 坐标处理

- mRNA 位点：`ATCG*.1` transcript ID 直接匹配 GTF 中 `ChrC` transcript 注释，
  按 exon 累计偏移换算到 `ChrC` genomic coordinate（负链反向）。
- rRNA（`Ch*`）与 tRNA（`tRNA-Pt-*`）位点表使用转录本内坐标；`plot.R` 内置
  `Ch*` / `tRNA-Pt-*` 到 `ATCG*` gene ID 的映射表，再根据 GTF strand 与 exon 结构换算。
- 输出的 `chloroplast_chr_overview_sites.tsv` 保留 `source_chr`、`gene_id`、`pos`、
  `chr_pos`、`ratio_pct` 等字段，便于检查映射结果。

## 当前结果

- ChrC feature：88 mRNA、8 rRNA、37 tRNA
- 映射到 ChrC 的 WT leaf 位点：4 rRNA、63 tRNA、18 mRNA
- 1 个 tRNA 位点无法映射：`tRNA-Pt-Ser-TGA-1:74`（超出该基因 exon 范围），
  与论文版一致，写入 `chloroplast_chr_overview_unmapped_sites.tsv` 并给出 warning

## 输出

- `chloroplast_chr_overview.pdf/png/svg`：环形总览图
- `chloroplast_chr_overview_annotated.pdf/png/svg`：同一张图 + mRNA 基因 ID 标注
- `chloroplast_chr_overview_features.tsv`：用于作图的 ChrC feature
- `chloroplast_chr_overview_sites.tsv`：映射到 ChrC 的 WT leaf 位点
- `chloroplast_chr_overview_unmapped_sites.tsv`：未能映射的位点
- `chloroplast_chr_overview_feature_summary.tsv`：feature 数量汇总
- `chloroplast_chr_overview_site_summary.tsv`：位点数与平均修饰比例汇总

图例：tRNA / mRNA / rRNA 三色同时用于外圈基因与内圈位点；径向 bar 长度正比于修饰比例，
黑色外点标记位点位置。

## 运行

```bash
bash data.sh
```
