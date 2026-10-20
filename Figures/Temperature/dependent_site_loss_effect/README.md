# dependent_site_loss_effect — trua1b dependent Ψ 位点的修饰下降效应散点图

## 目的

在 21℃ 和 27℃ 下**分别**看 TRUA1B 缺失造成的 Ψ 位点修饰下降：X 轴 = 突变体修饰
比例，Y 轴 = WT 修饰比例，灰点是该温度下全部 WT 阳性位点，红点是该温度下判出的
dependent 位点。红点落在对角线下方（WT > Mut）即为 TRUA1B 依赖的修饰下降。

移植自论文版 `archive/20260902/paper/TRUA1B_thermomorphogenesis/dependent_site_loss_effect/`。
三处不同：原图只有 mRNA，这里扩到 **mRNA / rRNA / tRNA 三类**（共 6 张）；数据源从旧
`data.sh` 现拼换成 2026-09-25 换血后的三重复 methylKit 口径；绘图样式改走 `Figures/shared`。

## 目录结构

```
Figures/Temperature/dependent_site_loss_effect/
├── data_prep.sh                              6 张表的构建逻辑（只产数据，不画图）
├── plot.R                                    逐张出图（pdf/png/svg）
└── <rna>_<temp>_dependent_site_loss_effect.tsv / .pdf / .png / .svg
```

## 数据来源

| RNA 类型 | 上游目录 |
|---|---|
| mRNA | `Temperature/STAT/dependent_methylKit/mRNA/<temp>/trua1b/` |
| rRNA | `Temperature/STAT/dependent_methylKit/rRNA/<temp>/trua1b/` |
| tRNA | `Temperature_tRNA/STAT/dependent_methylKit/ncRNA/<temp>/trua1b/` |

株系为 `trua1b_1` / `trua1b_3`（**trua1b_2 整线弃用**），每株系每温度 3 个生物学重复。

每个组合用到三张表：

- `trua1b_1_<temp>.dependent-sites.tsv` — 株系 1 的逐位点表，含**全部** WT 候选位点
- `trua1b_3_<temp>.dependent-sites.tsv` — 株系 3 同上
- `trua1b-<temp>.dependent.tsv` — 两株系 AND（交集）后的最终 dependent 位点

### 判定口径

见 `Temperature/STAT/dependent_methylKit/run.sh`（tRNA 侧为 `Temperature_tRNA` 同源脚本）：
methylKit `F` 检验、**BH q < 0.05**（方向为正，meth.diff > 0）、效应量
`Delta_Ratio(WT-Mut) >= 0.05`、Mock 平均深度 >= 10，且**两个株系同时成立**（AND）
才判为 dependent。

### 画图表口径

- 背景 = 该温度下的全部 WT 阳性位点（上游 fisher 的 `WT_<temp>` 阳性集）。
- `WT_Ratio` 取 per-line 表；`Mut_Ratio` = **两株系 Mut_Ratio 的均值**。已逐位点实测：
  该均值与 `trua1b-<temp>.dependent.tsv` 自带的 `Mut_Ratio` 完全一致（`data_prep.sh`
  因此可以用同一个公式统一处理全部位点，不必对 dependent 位点另走特例）。
- `site_status` = `dependent`（落在最终交集表内）/ `wt_only`（其余）。

### tRNA 侧的 ^tRNA- 过滤

`Temperature_tRNA` 的参考是 `tRNA-unique_sn_snoRNA.fa`，混有 snRNA / snoRNA，
**必须按 `^tRNA-` 筛掉**（与 `Figures/dependent_methylKit/tRNA`、同目录姊妹图
`temperature_responsive_scatter` 同口径）。上游 ncRNA 全表 614 / 608 个候选位点，
过滤后剩 503 / 508。

### 位点数

与 `Temperature/STAT/dependent_methylKit/summary.tsv`、`Temperature_tRNA` 同表逐一对齐：

| RNA 类型 | 21℃ WT 位点 | 21℃ dependent | 27℃ WT 位点 | 27℃ dependent |
|---|---|---|---|---|
| mRNA | 196 | 0 | 194 | 0 |
| rRNA | 106 | 0 | 107 | 0 |
| tRNA | 503 | 0 | 508 | 0 |

**六张图全灰**：2026-10-20 启用 BH q 后，三类 RNA 在两个温度下都没有判出
dependent 位点（旧版原始 p 口径下 mRNA 21℃ 有 5 个、tRNA 有 3 个）。按需求照出以
保持三类齐全，`plot.R` 遇到零 dependent 会打一条日志，免得日后误判成数据没接上。

## 运行

```bash
bash data_prep.sh     # 只产 6 张 TSV
Rscript plot.R        # 只画 6 张图
```

## 产物

6 套 `<rna>_<temp>_dependent_site_loss_effect.{pdf,png,svg}`，每套 6.8 × 6.8 in、
600 dpi，边框面板（`theme_base`），图例内嵌右下角，左上角标注温度。
