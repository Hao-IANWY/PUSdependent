# chloroplast_RNA — 叶绿体 Ψ landscape 与 TRUA4 依赖位点图组

## 目的

把研究范围从核编码转录本扩展到叶绿体 RNA：展示 WT 叶片叶绿体 rRNA / tRNA / mRNA 的
Ψ 位点数量与修饰比例分布，用叶绿体纯化 RNA 建库数据做正交验证，并展示 `trua4` 突变体
中叶绿体 tRNA 上的 TRUA4 依赖位点。

移植自论文版 `archive/20260902/paper/chloroplast_RNA/`，数据源换成当前流程，
绘图样式改走 `Figures/shared/plot_style.R` + `plot_tokens.R`（统一 Arial、无 showtext，
`save_plot()` 一次输出 PDF/PNG/SVG）。

## 目录结构

```
Figures/chloroplast_RNA/
├── chloroplast_site_overview/              WT leaf 叶绿体位点数量与比例分布
├── chloroplast_purified_rna_validation/    叶绿体纯化 RNA 位点数 + overlap Venn
├── trua4_dependent_site_overview/          WT vs trua4 叶绿体位点散点图
├── trua4_dependent_trna_heatmap/           trua4 dependent 叶绿体 tRNA 结构热图
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
| trua4 突变体 | `polyARNA/STAT/fisher/{rRNA,mRNA}/trua4_{1,2}-leaf_sites.Fisher.tsv`、`tRNA/STAT/fisher/ncRNA/trua4_{1,2}_Leaf_sites.Fisher.tsv` |
| trua4 dependent | `polyARNA/STAT/dependent_methylKit/{rRNA,mRNA}/leaf/trua4/`、`tRNA/STAT/dependent_methylKit/ncRNA/leaf/trua4/` |
| 叶绿体注释 | `Araport11_fixed.gtf`（`ChrC`） |

叶绿体位点的识别口径：

- rRNA：`chr == "Ch23S"`
- tRNA：`chr` 以 `tRNA-Pt-` 开头（排除同一 ncRNA 参考里的核/线粒体 tRNA 与 sn/snoRNA）
- mRNA：`chr` 以 `ATCG` 开头（叶绿体编码基因）

注：`Chloroplast/` 是指向测序流程工作目录的软链接，未入库；图中用到的派生表已一并提交。

## 与论文版实现的三处差异

1. **数据源更新**：论文版使用的 `temp_polyA` / `temp_totalRNA` / `chloroplast` 目录，
   现分别对应 `polyARNA` / `tRNA` / `Chloroplast`；位点表由当前流程重新生成。
2. **dependent 口径**：论文版用旧 Fisher/Welch dependent 表；这里统一用
   `STAT/dependent_methylKit`（methylKit 口径）。当前 per-line 表包含全部候选位点，
   绘图前按 `is_dependent == TRUE` 过滤。
3. **fix_pos 现算**：当前 dependent-sites 表不再带 `fix_pos`，热图脚本按
   `All-tRNA_rep-structure.csv` 的区段顺序与固定长度，把 native pos 重标定到
   canonical 0–76（与 `Figures/dependent_methylKit/tRNA/build_genotype_dependent_tables.sh` 一致）。

## 当前数据规模

| 图 | rRNA | tRNA | mRNA |
|---|---|---|---|
| WT leaf 位点（site overview / 散点图背景） | 4 | 64 | 18 |
| 叶绿体纯化 RNA 位点 | 3 | 59 | 15 |
| 两者 overlap | 3 | 55 | 11 |
| trua4 dependent（任一线） | 2 | 7 | 0 |

tRNA dependent 的 7 个位点中 6 个为两株系共享，1 个仅 `trua4_2`；热图最终落在
7 个 Pt tRNA 代表的反密码子茎（canonical 38–40）。

## 运行顺序

```bash
bash Figures/chloroplast_RNA/chloroplast_site_overview/data.sh
bash Figures/chloroplast_RNA/chloroplast_purified_rna_validation/data.sh  # 读取上一目录的产出表
bash Figures/chloroplast_RNA/trua4_dependent_site_overview/data.sh
bash Figures/chloroplast_RNA/trua4_dependent_trna_heatmap/data.sh
bash Figures/chloroplast_RNA/chloroplast_chr_overview/data.sh
```

`chloroplast_purified_rna_validation` 以 `chloroplast_site_overview` 的
`chloroplast_site_overview_density_sites.tsv` 作为 WT leaf 对照集，必须先跑前者。
其余目录彼此独立。
