# overlap — 27℃ 温度上调位点 vs 27℃ trua1b dependent 位点的 overlap Venn

## 目的

看 27℃ 下**修饰升高**的位点，和 27℃ 下**需要 TRUA1B** 的位点，是不是同一批位点。
左圈 = 27℃ 上调位点（WT_27 vs WT_21 双侧检验判为 up_in_WT_27），右圈 = 27℃ 的
trua1b dependent 位点，中间是交集计数。

移植自论文版 `archive/20260902/paper/TRUA1B_thermomorphogenesis/temperature_responsive_dependent_overlap/`
的 **`_27only` 2-set 变体**（原图只有 mRNA 一张，且是 2-set 与 3-set 一起出的），
这里只出 2-set 并按 RNA 类型拆成 mRNA / rRNA / tRNA 三张。

## 数据来源

| 集合 | mRNA / rRNA | tRNA |
|---|---|---|
| 27℃ 上调 | `Temperature/STAT/Temperature_response/<type>/WT_27_vs_WT_21.response.tsv`（`direction == up_in_WT_27`） | `Temperature_tRNA/STAT/Temperature_response/ncRNA/` 同名表 |
| 27℃ dependent | `Temperature/STAT/dependent_methylKit/<type>/27/trua1b/trua1b-27.dependent.tsv` | `Temperature_tRNA/STAT/dependent_methylKit/ncRNA/27/trua1b/` 同名表 |

位点主键 = `chr + pos + strand`；两组取并集后逐位点打标记，`combination` 记
`up_only` / `dependent_only` / `both`。

tRNA 侧两组都按 `^tRNA-` 过滤 snRNA / snoRNA（同 `Figures/Temperature` 其余图组口径）。
dependent 走 methylKit 口径，与同目录姊妹图 `dependent_site_loss_effect` 一致。

## 结果：三组交集全为 0

| RNA 类型 | 27℃ 上调 | 27℃ dependent | **交集** |
|---|---|---|---|
| mRNA | 6 | 0 | **0** |
| rRNA | 13 | 0 | **0** |
| tRNA | 0 | 0 | **0** |

2026-10-20 启用 BH q 后，mRNA 的 27℃ 上调位点从 10 降到 6，tRNA 从 3 降到 0
（Temperature_tRNA 的温度响应也为 0）；三类的右圈全为空，tRNA 左圈也为空。

**注意这是个阴性结果**，且与论文版不同：论文版同一分析（mRNA、2 重复、Fisher 口径）
是 26 / 8 / 1，交集有 8 个位点。新版 3 重复数据里 mRNA 在 27℃ 一个 dependent 位点都
没判出来，交集必然为 0。

换口径也不改变结论：`Temperature/STAT/dependent`（纯效应量版）在 27℃ 下 mRNA 有 4 个、
tRNA 有 4 个 dependent 位点，但**与上调集合仍然零重合**；放宽到 `pos ± 5bp` 容差同样
零重合。即：在本套数据里，「27℃ 修饰升高」和「TRUA1B 依赖」是两批互不相交的位点。

三类 RNA 在 BH 口径下 27℃ 都没有 dependent 位点，右圈全为空；按需求照出保持
三类齐全。

## 为什么手绘圆

`plot.R` 用 ggplot 的 `geom_polygon` 手画两个圆，而不是 `VennDiagram` 包
（仓库里 `Figures/dependent_methylKit/rRNA/sites_count/overlap` 用的是后者）。
手绘的好处正是**空集合也能照常出图**——`VennDiagram` 遇到空集合会直接失败，
而本分析里 mRNA / rRNA 的右圈恰恰是空的。圆的布局固定，不随集合大小或交集缩放。

## 运行

```bash
bash data_prep.sh     # 只产 3 张 overlap TSV
Rscript plot.R        # 只画 3 张 Venn
```

## 产物

- `<rna>_27up_vs_trua1b_dependent_overlap.tsv` — 逐位点标记表
- `<rna>_27up_vs_trua1b_dependent_venn.{pdf,png,svg}` — 9 × 7 in、600 dpi
