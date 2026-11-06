# temperature_responsive_scatter — WT 21/27 度温度响应 Ψ 位点散点图（tRNA）

## 目的

把 WT 在 21℃ 与 27℃ 两个温度下的 tRNA Ψ 位点摊在同一张散点图上：X 轴 = 21℃ 修饰
比例，Y 轴 = 27℃ 修饰比例，对角线为 y=x。落在对角线右下方的是 27℃ 修饰下降的位点
（蓝），左上方是上升的位点（红），其余为灰。

mRNA / rRNA 版在姊妹目录
`Figures/dependent_betabinomial/Temperature_mRNA/temperature_responsive_scatter/`。

## 目录结构

```
Figures/dependent_betabinomial/Temperature_tRNA/temperature_responsive_scatter/
├── data_prep.sh                              数据表构建逻辑（只产数据，不画图）
├── plot.R                                    出图（pdf/png/svg）
└── tRNA_temperature_responsive_scatter.tsv / .pdf / .png / .svg
```

图例固定三行顺序与配色（`drop = FALSE`），即使 down-regulated 位点为空图例仍齐全。

## 数据来源

唯一上游：
`Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.all-sites.tsv`

`all-sites.tsv` 就是**候选位点全表**（两个温度 fisher 阳性位点的并集，逐位点回
adjustGap 重算两个温度的修饰比与 Mock 平均深度），因此本目录不回原始
fisher/adjustGap 重算，直接取库内产物，只在 `data_prep.sh` 里做列名规范化与行筛选。

### 判定口径

见 `Temperature_tRNA/STAT/Temperature_response_betabinomial/run.sh`：

- 候选 = `WT_21_sites.Fisher.tsv` 与 `WT_27_sites.Fisher.tsv` 中 Pass==1 位点的并集；
  只在单一温度阳性的位点被刻意保留——那正是"另一温度修饰丢失"的强候选。
- beta-binomial GLMM 全局 theta + 逐位点 LRT、**BH q < 0.05**（双侧 FDR）、
  `|Delta_Ratio(WT_27-WT_21)| >= 0.05`、两个温度组每个 Mock/Treat 重复的 depth 都
  `>= 10`，才判为 temperature-responsive，方向记在 `direction`
  （`up_in_WT_27` / `down_in_WT_27`）。
- `data_prep.sh` 把 `direction` 空值改写为 `non-responsive`，另两值改写为
  `up-regulated` / `down-regulated`（= 27℃ 上调 / 下调）。

位点数与 `Temperature_tRNA/STAT/Temperature_response_betabinomial/summary.tsv`
逐一对齐（按 `^tRNA-` 过滤后口径）：

| RNA 类型 | 候选位点 | 响应位点 | 27℃ 上调 | 27℃ 下调 |
|---|---|---|---|---|
| tRNA | 541 | 0 | 0 | 0 |

输入为 2026-11-06 六重复重建（batch1 剔除，WT_21/WT_27 各 rep1~6 两批次口径）
的全量重跑结果。上游 ncRNA 全表 695 个候选（`summary.tsv`：695 candidate /
622 tested / 1 responsive，1 up / 0 down）；唯一响应位点是
`AT5G00750-sno:89`（snoRNA），按 `^tRNA-` 过滤后 tRNA 保留 541 个候选
（154 个 sn/snoRNA 候选被筛掉）、**0 个响应位点**——散点图为全灰。

### tRNA 侧的 ^tRNA- 过滤

`Temperature_tRNA` 的参考是 `tRNA-unique_sn_snoRNA.fa`，混有 snRNA / snoRNA，
**必须按 `^tRNA-` 筛掉**（与 `Figures/dependent_betabinomial/tRNA` 同口径）。
当前上游 ncRNA 全表 695 个候选，过滤后剩 541 个（154 个 sn/snoRNA 候选被筛掉）。

## 运行

```bash
bash data_prep.sh     # 只产 TSV
Rscript plot.R        # 只画图
```

## 产物

`tRNA_temperature_responsive_scatter.{pdf,png,svg}`，6.8 × 6.8 in、600 dpi，
边框面板（`theme_base`），图例内嵌右下角。
