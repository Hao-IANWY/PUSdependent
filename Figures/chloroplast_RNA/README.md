# chloroplast_RNA — 叶绿体 Ψ landscape 图组

## 目的

把研究范围从核编码转录本扩展到叶绿体 RNA：展示 WT 叶片叶绿体 rRNA / tRNA / mRNA 的
Ψ 位点数量与修饰比例分布，并用叶绿体纯化 RNA 建库数据做正交验证。

移植自论文版 `archive/20260902/paper/chloroplast_RNA/`，数据源换成当前流程，
绘图样式改走 `Figures/shared/plot_style.R` + `plot_tokens.R`（统一 Arial、无 showtext，
`save_plot()` 一次输出 PDF/PNG/SVG）。

注：trua4 dependent 散点图与 Pt tRNA 结构热图两个子目录原以
`STAT/dependent_methylKit`（methylKit 口径）为数据源，methylKit 版结果
2026-10-30 整体移除后一并下线；如需恢复可基于
`STAT/dependent_betabinomial` 的 per-line 表（`is_dependent == TRUE`）重建。

## 目录结构

```
Figures/chloroplast_RNA/
├── chloroplast_site_overview/              WT leaf 叶绿体位点数量与比例分布
├── chloroplast_purified_rna_validation/    叶绿体纯化 RNA 位点数 + overlap Venn
├── chloroplast_chr_overview/               ChrC 环形染色体总览（附 mRNA 基因标注版）
└── README.md
```

## 数据源

| 数据 | 路径 |
|---|---|
| WT leaf rRNA | `polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv` |
| WT leaf mRNA | `polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv` |
| WT leaf tRNA | `tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv` |
| 叶绿体纯化 RNA | `Chloroplast/STAT/fisher/{ncRNA,mRNA}/Chloroplast_sites.positive.tsv` |
| 叶绿体注释 | `Araport11_fixed.gtf`（`ChrC`） |

叶绿体位点的识别口径：

- rRNA：`chr == "Ch23S"`
- tRNA：`chr` 以 `tRNA-Pt-` 开头（排除同一 ncRNA 参考里的核/线粒体 tRNA 与 sn/snoRNA）
- mRNA：`chr` 以 `ATCG` 开头（叶绿体编码基因）

注：`Chloroplast/` 是指向测序流程工作目录的软链接，未入库；图中用到的派生表已一并提交。

## 与论文版实现的差异

1. **数据源更新**：论文版使用的 `temp_polyA` / `temp_totalRNA` / `chloroplast` 目录，
   现分别对应 `polyARNA` / `tRNA` / `Chloroplast`；位点表由当前流程重新生成。

## 当前数据规模

| 图 | rRNA | tRNA | mRNA |
|---|---|---|---|
| WT leaf 位点（site overview） | 4 | 64 | 18 |
| 叶绿体纯化 RNA 位点 | 3 | 59 | 15 |
| 两者 overlap | 3 | 55 | 11 |

## 运行顺序

```bash
bash Figures/chloroplast_RNA/chloroplast_site_overview/data.sh
bash Figures/chloroplast_RNA/chloroplast_purified_rna_validation/data.sh  # 读取上一目录的产出表
bash Figures/chloroplast_RNA/chloroplast_chr_overview/data.sh
```

`chloroplast_purified_rna_validation` 以 `chloroplast_site_overview` 的
`chloroplast_site_overview_density_sites.tsv` 作为 WT leaf 对照集，必须先跑前者。
其余目录彼此独立。
