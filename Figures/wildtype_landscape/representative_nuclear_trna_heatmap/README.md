# representative_nuclear_trna_heatmap

## 目的
Panel E「tRNA landscape」：展示代表性核 tRNA 的 positive-site ratio 沿 **canonical tRNA 二级结构布局**的分布，三个组织各一张（Root 为主图，Leaf / Flower 同脚本一并产出）。

- **行** = 46 个核 tRNA 家族（代表成员 = `family_top` 里 mean_CPM 最高的那条 tRNA）
- **列** = canonical 0–76 位（按 p1 / Accstem / Dstem / Dloop / Anticodon stem-loop / Vregion / Tstem / Tloop / CCA 分块，底部有区段色带、位号与区段名）
- **值** = 该家族在该结构位上**所有拷贝位点**的 depth 加权平均 ratio；无位点的单元格为浅灰（`#FCFCFC`）
- **配色** = `#FFF6EB → #C44E00`（Ratio 0–1）

## 代表 tRNA 与 Sprinzl 坐标的来源

- **代表 tRNA**：`tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`。该表每个 anticodon family
  （分室 + 氨基酸 + 反密码子）一行，`feature_id` 是该家族里 `mean_CPM`（Mock 全部文库均值）最高的成员。
  本图取其中的 **Nu 家族 46 个**，行序按 `family` 字典序；代表成员记在映射表的 `Representative_ID` 列。
- **Sprinzl 坐标**：直接取自 `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv`
  的 `Sprinzl_pos` 列（由 `BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py`
  基于 `plantRNA_database_unique.tsv` 的 220 条 complete 注释生成）。
  **本目录不做任何坐标转换**，`plot.R` 只把 canonical 位号查表成展示列。
- 每个拷贝用**自己**的区段注释，因此不再有“拷贝长度与家族代表不符就剔除”的逻辑。
- 只画核 tRNA；Mt / Pt 的 tRNA 与 sn/snoRNA 行不计入（汇总表里计为 `excluded_nonNu`）。

## 运行流程
```bash
bash data.sh          # = 折叠(awk) + 画图(Rscript plot.R)
```

配图为论文版实现
`archive/20260902/paper/WT_sites_overview/representative_nuclear_trna_heatmap/plot_representative_nuclear_trna_heatmap.R`
的移植；投影改为 `data.sh` 直接读取已注释的 `Sprinzl_pos`，绘图与样式走
`Figures/shared/plot_style.R` + `plot_tokens.R`（`save_plot()` 一次出 PDF/PNG/SVG，字体统一 Arial、无 showtext）。

## 主要输入
- `../../../tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv`
  （列：`chr, pos, Sprinzl_pos, strand, ratio, depth`）
- `../../../tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`（家族代表，表达量口径）

## 数据接口

`data.sh` 生成三个映射表，作为 `plot.R` 的唯一输入：

`representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv`，列为
`tissue, Rep_ID, Representative_ID, chr, pos, strand, ratio, depth, canonical_pos`

- `Rep_ID` = 家族名（如 `Nu-Arg-TCG`），同一家族的所有拷贝位点都并入该行；
- `Representative_ID` = `family_top` 里该家族表达量最高的拷贝（如 `tRNA-Nu-Arg-TCG-2`）；
- `canonical_pos` = 已注释的 `Sprinzl_pos`（0–76）；
- 同一 (家族, canonical 位) 收到多个原始位点时按 depth 加权合并。

当前规模（`trna_landscape_summary.tsv`）：

| tissue | 读入行 | 非核/非 tRNA | Sprinzl 缺失 | 映射位点 | 覆盖家族 | 有值单元格 |
|---|---|---|---|---|---|---|
| Root | 439 | 84 | 0 | 355 | 46 | 126 |
| Leaf | 567 | 167 | 0 | 400 | 46 | 134 |
| Flower | 550 | 154 | 0 | 396 | 46 | 132 |

## 与旧版（区段表投影）的差异

旧版用 `Nu-tRNA_rep-structure.csv` / `Nu-tRNA_rep_structure_regions.tsv`（论文版策展的 46 条代表）
在 `data.sh` 里做重标定。改成当前方案后：

- 同一位点对拍：Root 旧版会纳入的 350 个位点中 **348 个一致**；Leaf 395 中 393 一致；Flower 391 中 389 一致。
- 不一致的 2 个位点（每组织相同）：`tRNA-Nu-Gly-GCC-2` / `-5` 的 pos 13，旧版给 13、新版给 **14**。
  这两个拷贝的 D 臂区段划分与家族代表不同（3+9+3 vs 4+8+4），按拷贝自己的注释应为 14。
- 旧版剔除的 3 个长度异常拷贝（`tRNA-Nu-Asn-GTT-5`、`tRNA-Nu-Leu-TAA-3`、`tRNA-Nu-Ser-AGA-5`）
  共 5 个位点/组织，新版按拷贝注释正常纳入。
- 旧版依赖的 `Nu-tRNA_rep-structure.csv` / `Nu-tRNA_rep_structure_regions.tsv` 已删除，
  该目录不再保留任何区段注释表（统一由 `plantRNA_database_unique.tsv` 提供）。

## 输出
- `representative_nuclear_trna_heatmap_{root,leaf,flower}.pdf/png/svg`（各 18 in 宽 × 11.08 in 高，600 dpi）
- `representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv`
- `trna_landscape_summary.tsv`

## 验证
- 与旧映射逐位点对拍：一致率 348/350、393/395、389/391，差异全部落在上述可解释的拷贝上。
- 字体：PDF 内嵌 Arial-BoldMT 子集，SVG 全部 `font-family: "Arial"`，无 DejaVu。
- 生物学 sanity：T loop 位（canonical 54–60）出现纵向高 ratio 带，与论文版图一致。

## 已知无害提示
`ComplexHeatmap::Heatmap()` 构造对象时会在"无设备"状态下测量文本宽度，R 临时打开的默认 pdf 设备
找不到 Arial，会抛一批 `PostScript字体数据库里找不到'Arial'` 警告。测量值与 cairo 设备下一致，
且"构造一次、三种格式共用"保证 PDF/PNG/SVG 布局完全一致，故 `plot.R` 只定向屏蔽这一条警告，其他警告照常抛出。
