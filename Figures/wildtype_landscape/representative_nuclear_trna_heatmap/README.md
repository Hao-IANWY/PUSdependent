# representative_nuclear_trna_heatmap

## 目的
Panel E「tRNA landscape」：展示代表性核 tRNA 的 positive-site ratio 沿 **canonical tRNA 二级结构布局**的分布，三个组织各一张（Root 为主图，Leaf / Flower 同脚本一并产出）。

- **行** = 46 个核 tRNA 代表
- **列** = canonical 0–76 位（按 p1 / Accstem / Dstem / Dloop / Anticodon stem-loop / Vregion / Tstem / Tloop / CCA 分块，底部有区段色带、位号与区段名）
- **值** = 该代表在该结构位上所有拷贝位点的 **depth 加权平均 ratio**；无位点的单元格为浅灰（`#FCFCFC`）
- **配色** = `#FFF6EB → #C44E00`（Ratio 0–1）

## 代表性 tRNA 与结构布局的来源

- 行 = `tRNA/count/ncRNA/family_top` 的 **Nu 家族**（基因组+氨基酸+反密码子层，共 46 个），行序沿用 `Nu-tRNA_rep-structure.csv`。
- `Nu-tRNA_rep-structure.csv` 与 `Nu-tRNA_rep_structure_regions.tsv` 直接取自
  `archive/20260902/paper/WT_sites_overview/representative_nuclear_trna_heatmap/`（论文版实现的结构注释）。
- **坐标兼容性已实测**：archive 结构序列与 `tRNA/config.yaml` 所用参考
  `/home/HL/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa` 逐碱基一致（每个家族恰有一个拷贝完全匹配），
  结构表坐标跨度 = 序列长度 ⇒ 结构表坐标就是参考 fasta 的 1-based 坐标（**含 CCA 尾**）；
  三个组织的全部新 positive site 都落在结构跨度内（46/46 家族、零越界）。因此**无需重做 tRNAscan/区段注释**。
- 只画核 tRNA。Pt / Mt 不画：环境无 tRNAscan-SE / infernal，且质体/线粒体 tRNA 的环长分布与核基因不同，投影到同一 canonical 布局有扭曲风险。

## 运行流程
```bash
bash data.sh          # = 映射(awk) + 画图(Rscript plot.R)
```

配图为论文版实现
`archive/20260902/paper/WT_sites_overview/representative_nuclear_trna_heatmap/plot_representative_nuclear_trna_heatmap.R`
的移植：映射改为 `data.sh` 前置完成，绘图与样式改走 `Figures/shared/plot_style.R` + `plot_tokens.R`（`save_plot()` 一次出 PDF/PNG/SVG，字体统一 Arial、无 showtext）。

## 主要输入
- `../../../tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.tsv`
- `Nu-tRNA_rep-structure.csv`、`Nu-tRNA_rep_structure_regions.tsv`（见上）
- `/home/HL/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa`（仅用于识别长度与代表不符的拷贝）

## 数据接口
`data.sh` 生成三个映射表，作为 `plot.R` 的唯一输入：

`representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv`，列为
`tissue, Rep_ID, chr, pos, strand, ratio, depth, region, region_pos, fixed_pos, matrix_col_id, canonical_pos`

投影规则：位点落在区段 `[start_pos, end_pos]` 内 → `region_pos = pos - start_pos + 1`，再按区段固定长度重标定
`round((region_pos - 1) * (fixed_length - 1) / (region_length - 1)) + 1`（`round` 为 R 的银行家舍入，必须在 awk 里保持一致）；
`canonical_pos` 为 0–76 的展示位号。同一 canonical 位上收到多个原始位点时按 depth 加权合并。

当前规模（`trna_landscape_summary.tsv`）：

| tissue | 读入 tRNA-Nu 位点 | 排除长度不符拷贝 | 未命中区段 | 映射位点 | 覆盖代表 | 有值单元格 |
|---|---|---|---|---|---|---|
| Root | 370 | 5 | 0 | 365 | 46 | 128 |
| Leaf | 407 | 5 | 0 | 402 | 46 | 137 |
| Flower | 405 | 5 | 0 | 400 | 46 | 137 |

**排除说明**：3 个拷贝的参考序列长度与所属家族代表不一致（`tRNA-Nu-Asn-GTT-5`、`tRNA-Nu-Leu-TAA-3`、`tRNA-Nu-Ser-AGA-5`），
坐标会整体位移，故剔除；每组织 5 个位点（Asn-GTT-5 的 57；Leu-TAA-3 的 40、66；Ser-AGA-5 的 41、66）。

## 输出
- `representative_nuclear_trna_heatmap_{root,leaf,flower}.pdf/png/svg`（各 18 in 宽 × 11.08 in 高，600 dpi）
- `representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv`
- `trna_landscape_summary.tsv`

## 验证
- 独立复核：python 从原始 fisher 表独立重算投影与加权均值、base R 从映射表独立折叠，两个实现与产出逐格对拍
  （三组织键集合完全一致，最大数值差 3.33e-11）。
- 字体：PDF 内嵌 Arial-BoldMT 子集，SVG 全部 `font-family: "Arial"`，无 DejaVu。
- 生物学 sanity：T loop 位（canonical 54–60）出现纵向高 ratio 带，与论文版图一致。

## 已知无害提示
`ComplexHeatmap::Heatmap()` 构造对象时会在"无设备"状态下测量文本宽度，R 临时打开的默认 pdf 设备
找不到 Arial，会抛一批 `PostScript字体数据库里找不到'Arial'` 警告。测量值与 cairo 设备下一致，
且"构造一次、三种格式共用"保证 PDF/PNG/SVG 布局完全一致，故 `plot.R` 只定向屏蔽这一条警告，其他警告照常抛出。
