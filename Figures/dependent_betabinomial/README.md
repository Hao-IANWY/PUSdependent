# dependent_betabinomial — beta-binomial 口径的 dependent 位点图组

`Figures/dependent_methylKit` 的姊妹图组：脚本逐目录移植（绘图实现、样式、
目录结构与运行方式完全一致），数据源全部从 `STAT/dependent_methylKit` 换成
`STAT/dependent_betabinomial`（`call_dependent_betabinomial.R`，Jie Huang
beta-binomial GLMM 方法；见仓库 `dependent_jiehuang/`）。两版本的判定门槛
（BH q < 0.05、Δ ≥ 0.05、深度门、两株系 AND）一致，只有统计检验不同。

```
Figures/dependent_betabinomial/
├── tRNA/    tRNA 图组全量：sites_count(+overlap) / loss_effect / motif /
│            merge_heatmap / 7 个基因型热图+分布箱线图 / count_missed_dependent
├── mRNA/    polyARNA：sites_count(+overlap) / loss_effect / ratio_distribution /
│            region_distribution
└── rRNA/    polyARNA：sites_count(+overlap) / loss_effect / tissue_specificity
```

## 与 methylKit 图组的差异

- **overlap Venn 与 methylKit 图组同构**：对比「beta-binomial（本图组方法，
  红）vs 纯效应量版 dependent（蓝）」；beta-binomial 的位点是其完全子集
  （544⊂952 / 53⊂299 / 21⊂44），与 methylKit 版「完全子集」结论一致。
  beta-binomial vs methylKit 的逐位点对照见各 STAT 目录，不在这张图里。
- **count_missed_dependent.sh 适配了行表列位**：beta-binomial 行表的
  status/is_dependent 在 14/15 列（多出 beta_coef/beta_se/odds_ratio 三列），
  方向判据用 `beta_coef < 0`（≙ methylKit 的 `meth.diff > 0`）；深度中位数
  改从 adjustGap 按候选宇宙键回查（beta-binomial 无 methylkit.samples 中间件）。
- **未移植**：`schematic/`（methylKit 检验流程示意图，属 methylKit 论文版）、
  `mRNA/rRNA 的 trua1a_dependent` 深挖目录、`tRNA/trua3/root_depth`
  （methylKit Root 无功效的机制归因，beta-binomial 版无此现象）。

## 计数快照（2026-10-25）

| 图组 | beta-binomial | methylKit | 交集 | 备注 |
|---|---|---|---|---|
| tRNA（^tRNA- 过滤，21 组合） | 544 | 468 | 397 | Root 154 vs 44：methylKit line2 φ 爆炸（q 齐卡 0.055）被 beta-binomial 全局 theta 救回 |
| polyARNA mRNA | 53 | 81 | 47 | beta-binomial 的 BH 更保守 |
| polyARNA rRNA | 21 | 17 | 10 | |

各图组 data_prep.sh 可随时从对应 `STAT/dependent_betabinomial` 重算。
