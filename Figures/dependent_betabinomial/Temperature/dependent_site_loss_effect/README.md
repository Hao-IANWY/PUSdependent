# dependent_site_loss_effect — trua1b dependent Ψ 位点的修饰下降效应散点图

## 目的

在 21℃ 和 27℃ 下**分别**看 TRUA1B 缺失造成的 Ψ 位点修饰下降：X 轴 = 突变体修饰
比例，Y 轴 = WT 修饰比例，灰点是该温度下全部 WT 阳性位点，红点是该温度下判出的
dependent 位点。红点落在对角线下方（WT > Mut）即为 TRUA1B 依赖的修饰下降。

移植自论文版 `archive/20260902/paper/TRUA1B_thermomorphogenesis/dependent_site_loss_effect/`。
三处不同：原图只有 mRNA，这里扩到 **mRNA / rRNA / tRNA 三类**（共 6 张）；数据源从旧
`data.sh` 现拼换成三重复 beta-binomial LRT 口径（mRNA / rRNA 为 2026-09-25 换血版，
tRNA 为 2026-11-02 全量换血版）；绘图样式改走 `Figures/shared`。

## 目录结构

```
Figures/dependent_betabinomial/Temperature/dependent_site_loss_effect/
├── data_prep.sh                              6 张表的构建逻辑（只产数据，不画图）
├── plot.R                                    逐张出图（pdf/png/svg）
└── <rna>_<temp>_dependent_site_loss_effect.tsv / .pdf / .png / .svg
```

## 数据来源

| RNA 类型 | 上游目录 |
|---|---|
| mRNA | `Temperature/STAT/dependent_betabinomial/mRNA/<temp>/trua1b/` |
| rRNA | `Temperature/STAT/dependent_betabinomial/rRNA/<temp>/trua1b/` |
| tRNA | `Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/<temp>/trua1b/` |

株系为 `trua1b_1` / `trua1b_3`（**trua1b_2 整线弃用**），每株系每温度 3 个生物学重复。
mRNA / rRNA 输入自 2026-09-25 换血后未变；tRNA 输入为 2026-11-02 全量换血后的
`Temperature_tRNA/STAT/dependent_betabinomial` 重跑结果（2026-11-02 14:54 完成）。

每个组合用到三张表：

- `trua1b_1_<temp>.dependent-sites.tsv` — 株系 1 的逐位点表，含**全部** WT 候选位点
- `trua1b_3_<temp>.dependent-sites.tsv` — 株系 3 同上
- `trua1b-<temp>.dependent.tsv` — 两株系 AND（交集）后的最终 dependent 位点

### 判定口径

见 `Temperature/STAT/dependent_betabinomial/run.sh`（tRNA 侧为 `Temperature_tRNA` 同源脚本）：
beta-binomial LRT `F` 检验、**BH q < 0.05**（方向为正，meth.diff > 0）、效应量
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
**必须按 `^tRNA-` 筛掉**（与 `Figures/dependent_betabinomial/tRNA`、同目录姊妹图
`temperature_responsive_scatter` 同口径）。当前上游 ncRNA 全表 547 / 555 个候选位点
（21℃ / 27℃），过滤后剩 493 / 493。

### 位点数

与 `Temperature/STAT/dependent_betabinomial/summary.tsv`、`Temperature_tRNA` 同表逐一对齐：

| RNA 类型 | 21℃ WT 位点 | 21℃ dependent | 27℃ WT 位点 | 27℃ dependent |
|---|---|---|---|---|
| mRNA | 196 | 0 | 194 | 0 |
| rRNA | 106 | 0 | 107 | 0 |
| tRNA | 493 | 0 | 493 | 0 |

**六张图全灰**：2026-10-20 启用 BH q 后，三类 RNA 在两个温度下都没有判出
dependent 位点（旧版原始 p 口径下 mRNA 21℃ 有 5 个、tRNA 有 3 个）。2026-11-02
tRNA 全量换血重跑后仍为 0/0（WT 位点数由旧版 503 / 508 变为 493 / 493）。按需求
照出以保持三类齐全，`plot.R` 遇到零 dependent 会打一条日志，免得日后误判成数据
没接上。

## 运行

```bash
bash data_prep.sh     # 只产 6 张 TSV
Rscript plot.R        # 只画 6 张图
```

## 产物

6 套 `<rna>_<temp>_dependent_site_loss_effect.{pdf,png,svg}`，每套 6.8 × 6.8 in、
600 dpi，边框面板（`theme_base`），图例内嵌右下角，左上角标注温度。
