# dependent_site_loss_effect — trua1b dependent Ψ 位点的修饰下降效应散点图（tRNA）

## 目的

在 21℃ 和 27℃ 下**分别**看 TRUA1B 缺失造成的 tRNA Ψ 位点修饰下降：X 轴 = 突变体
修饰比例，Y 轴 = WT 修饰比例，灰点是该温度下全部 WT 阳性位点，红点是该温度下判出的
dependent 位点。红点落在对角线下方（WT > Mut）即为 TRUA1B 依赖的修饰下降。

mRNA / rRNA 版在姊妹目录
`Figures/dependent_betabinomial/Temperature_mRNA/dependent_site_loss_effect/`。

## 目录结构

```
Figures/dependent_betabinomial/Temperature_tRNA/dependent_site_loss_effect/
├── data_prep.sh                              2 张表的构建逻辑（只产数据，不画图）
├── plot.R                                    逐张出图（pdf/png/svg）
└── tRNA_<temp>_dependent_site_loss_effect.tsv / .pdf / .png / .svg
```

## 数据来源

唯一上游：`Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/<temp>/trua1b/`

株系为 `trua1b_1` / `trua1b_2`（原 trua1b_3 于 2026-11-04 正名为 trua1b_2，与诺禾旧批
35 系同系），每株系（含 WT）6 个生物学重复。输入为 2026-11-06 六重复重建
（batch1 剔除，rep1~6 两批次口径）后由 `Temperature_tRNA/STAT/dependent_betabinomial`
重跑产出的结果。

每个组合用到三张表：

- `trua1b_1_<temp>.dependent-sites.tsv` — 株系 1 的逐位点表，含**全部** WT 候选位点
- `trua1b_2_<temp>.dependent-sites.tsv` — 株系 2（原株系 3）同上
- `trua1b-<temp>.dependent.tsv` — 两株系 AND（交集）后的最终 dependent 位点

### 判定口径

见 `Temperature_tRNA/STAT/dependent_betabinomial/run.sh`：beta-binomial LRT 检验、
**BH q < 0.05**（方向为正，meth.diff > 0）、效应量 `Delta_Ratio(WT-Mut) >= 0.05`、
Mock 平均深度 >= 10，且**两个株系同时成立**（AND）才判为 dependent。

### 画图表口径

- 背景 = 该温度下的全部 WT 阳性位点（上游 fisher 的 `WT_<temp>` 阳性集）。
- `WT_Ratio` 取 per-line 表；`Mut_Ratio` = **两株系 Mut_Ratio 的均值**。已逐位点实测：
  该均值与 `trua1b-<temp>.dependent.tsv` 自带的 `Mut_Ratio` 完全一致（`data_prep.sh`
  因此可以用同一个公式统一处理全部位点，不必对 dependent 位点另走特例）。
- `site_status` = `dependent`（落在最终交集表内）/ `wt_only`（其余）。

### tRNA 侧的 ^tRNA- 过滤

`Temperature_tRNA` 的参考是 `tRNA-unique_sn_snoRNA.fa`，混有 snRNA / snoRNA，
**必须按 `^tRNA-` 筛掉**（与 `Figures/dependent_betabinomial/tRNA`、姊妹目录
`temperature_responsive_scatter` 同口径）。当前上游 21℃ / 27℃ 全表 662 / 676 个
候选位点，过滤后剩 528 / 535。

### 位点数

与 `Temperature_tRNA/STAT/dependent_betabinomial/summary.tsv` 逐一对齐：

| RNA 类型 | 21℃ WT 位点 | 21℃ dependent | 27℃ WT 位点 | 27℃ dependent |
|---|---|---|---|---|
| tRNA | 528 | 0 | 535 | 0 |

**两张图全灰**：六重复新口径下两个温度都没有判出 dependent 位点。按需求照出以
保持口径完整，`plot.R` 遇到零 dependent 会打一条日志，免得日后误判成数据没接上。

## 运行

```bash
bash data_prep.sh     # 只产 2 张 TSV
Rscript plot.R        # 只画 2 张图
```

## 产物

2 套 `tRNA_<temp>_dependent_site_loss_effect.{pdf,png,svg}`，每套 6.8 × 6.8 in、
600 dpi，边框面板（`theme_base`），图例内嵌右下角，左上角标注温度。
