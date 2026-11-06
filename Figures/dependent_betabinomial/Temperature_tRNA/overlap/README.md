# overlap — 27℃ 温度上调位点 vs 27℃ trua1b dependent 位点的 overlap Venn（tRNA）

## 目的

看 tRNA 在 27℃ 下**修饰升高**的位点，和 27℃ 下**需要 TRUA1B** 的位点，是不是同一批
位点。左圈 = 27℃ 上调位点（WT_27 vs WT_21 双侧检验判为 up_in_WT_27），右圈 = 27℃
的 trua1b dependent 位点，中间是交集计数。

mRNA / rRNA 版在姊妹目录 `Figures/dependent_betabinomial/Temperature_mRNA/overlap/`。

## 数据来源

| 集合 | 上游 |
|---|---|
| 27℃ 上调 | `Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.response.tsv`（`direction == up_in_WT_27`） |
| 27℃ dependent | `Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/27/trua1b/trua1b-27.dependent.tsv` |

位点主键 = `chr + pos + strand`；两组取并集后逐位点打标记，`combination` 记
`up_only` / `dependent_only` / `both`。

两组都按 `^tRNA-` 过滤 snRNA / snoRNA（与 `Figures/dependent_betabinomial/tRNA`
同口径）。dependent 走 beta-binomial LRT 口径，与同目录姊妹图
`dependent_site_loss_effect` 一致。tRNA 输入为 2026-11-06 六重复重建（batch1 剔除，
rep1~6 两批次口径）后由 `Temperature_tRNA/STAT` 两个同名流程重跑产出的结果。

## 结果：三组计数全为 0

| RNA 类型 | 27℃ 上调 | 27℃ dependent | **交集** |
|---|---|---|---|
| tRNA | 0 | 0 | **0** |

2026-11-06 六重复重建后，tRNA 在 27℃ 没有判出任何温度上调位点（上游 ncRNA 全表
唯一响应位点是 `AT5G00750-sno:89`，snoRNA，被 `^tRNA-` 过滤筛掉），也没有
dependent 位点——两个圈都是空的，交集必然为 0。

## 为什么手绘圆

`plot.R` 用 ggplot 的 `geom_polygon` 手画两个圆，而不是 `VennDiagram` 包
（仓库里 `Figures/dependent_betabinomial/rRNA/sites_count/overlap` 用的是后者）。
手绘的好处正是**空集合也能照常出图**——`VennDiagram` 遇到空集合会直接失败，
而本分析里两个圈都是空的（2026-11-06 六重复口径）。圆的布局固定，不随
集合大小或交集缩放。

## 运行

```bash
bash data_prep.sh     # 只产 overlap TSV
Rscript plot.R        # 只画 Venn
```

## 产物

- `tRNA_27up_vs_trua1b_dependent_overlap.tsv` — 逐位点标记表
- `tRNA_27up_vs_trua1b_dependent_venn.{pdf,png,svg}` — 9 × 7 in、600 dpi
