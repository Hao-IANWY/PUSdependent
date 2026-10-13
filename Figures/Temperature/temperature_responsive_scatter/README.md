# temperature_responsive_scatter — WT 21/27 度温度响应 Ψ 位点散点图（mRNA / rRNA / tRNA）

## 目的

把 WT 在 21℃ 与 27℃ 两个温度下的 Ψ 位点摊在同一张散点图上：X 轴 = 21℃ 修饰比例，
Y 轴 = 27℃ 修饰比例，对角线为 y=x。落在对角线右下方的是 27℃ 修饰下降的位点
（蓝），左上方是上升的位点（红），其余为灰。

移植自论文版 `archive/20260902/paper/TRUA1B_thermomorphogenesis/temperature_responsive_sites/`，
三处不同：原图只有 mRNA，这里按 RNA 类型拆成 **mRNA / rRNA / tRNA 三张**；数据源从旧
Fisher 口径换成 2026-09-25 换血后的三重复 methylKit 口径；绘图样式改走 `Figures/shared`。

## 目录结构

```
Figures/Temperature/temperature_responsive_scatter/
├── data_prep.sh                              三张表的构建逻辑（只产数据，不画图）
├── plot.R                                    逐张出图（pdf/png/svg）
└── <rna>_temperature_responsive_scatter.tsv / .pdf / .png / .svg
```

三张图共用同一套图例顺序与配色（`drop = FALSE`），tRNA 没有 down-regulated 位点时
图例仍保留三行，三张图可直接并排比较。

## 数据来源

两个 STAT 流程，产物列名完全一致，故 `data_prep.sh` 一份逻辑跑三遍：

| RNA 类型 | 唯一上游 |
|---|---|
| mRNA | `Temperature/STAT/Temperature_response/mRNA/WT_27_vs_WT_21.all-sites.tsv` |
| rRNA | `Temperature/STAT/Temperature_response/rRNA/WT_27_vs_WT_21.all-sites.tsv` |
| tRNA | `Temperature_tRNA/STAT/Temperature_response/ncRNA/WT_27_vs_WT_21.all-sites.tsv` |

`all-sites.tsv` 就是**候选位点全表**（两个温度 fisher 阳性位点的并集，逐位点回
adjustGap 重算两个温度的修饰比与 Mock 平均深度），与 archive 版由 `data.sh` 现拼的
`temperature_responsive_scatter.tsv` 语义等价，因此本目录不回原始 fisher/adjustGap
重算，直接取库内产物，只在 `data_prep.sh` 里做列名规范化与行列筛选。

### 判定口径

见 `Temperature/STAT/Temperature_response/run.sh`（tRNA 侧为 `Temperature_tRNA` 同源脚本）：

- 候选 = `WT_21_sites.Fisher.tsv` 与 `WT_27_sites.Fisher.tsv` 中 Pass==1 位点的并集；
  只在单一温度阳性的位点被刻意保留——那正是"另一温度修饰丢失"的强候选。
- methylKit `F` 检验**双侧** p < 0.05、`|Delta_Ratio(WT_27-WT_21)| >= 0.05`、
  两个温度的 Mock 平均深度都 `>= 10`，才判为 temperature-responsive，方向记在
  `direction`（`up_in_WT_27` / `down_in_WT_27`）。
- `data_prep.sh` 把 `direction` 空值改写为 `non-responsive`，另两值改写为
  `up-regulated` / `down-regulated`（= 27℃ 上调 / 下调）。

位点数与各 `STAT/Temperature_response/summary.tsv` 逐一对齐：

| RNA 类型 | 候选位点 | 响应位点 | 27℃ 上调 | 27℃ 下调 |
|---|---|---|---|---|
| mRNA | 267 | 51 | 10 | 41 |
| rRNA | 108 | 25 | 13 | 12 |
| tRNA | 517 | 3 | 3 | 0 |

### tRNA 侧的 ^tRNA- 过滤

`Temperature_tRNA` 的参考是 `tRNA-unique_sn_snoRNA.fa`，混有 snRNA / snoRNA，
**必须按 `^tRNA-` 筛掉**（与 `Figures/dependent_methylKit/tRNA` 同口径）。上游
ncRNA 全表 638 个候选、6 个响应位点，其中 3 个 snoRNA 被过滤后剩 517 / 3。

## 运行

```bash
bash data_prep.sh     # 只产三张 TSV
Rscript plot.R        # 只画图
```

## 产物

`mRNA_` / `rRNA_` / `tRNA_` 三套 `temperature_responsive_scatter.{pdf,png,svg}`，
每套 6.8 × 6.8 in、600 dpi，边框面板（`theme_base`），图例内嵌右下角。
