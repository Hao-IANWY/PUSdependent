# tissue_ratio/tRNA — Flower 修饰比例偏低的 tRNA 位点

## 目的
在 WT Root、Leaf、Flower 三个组织之间比较 tRNA（ncRNA 参考中的 `^tRNA-` 位点）的修饰比例，
用 beta-binomial 统计检验找出 **Flower 中修饰比例显著偏低** 的位点，并沿代表性核 tRNA
家族 × canonical 0–76 结构布局展示（与 `representative_nuclear_trna_heatmap` 同构）。

## 方法（BID-transcript.script/downstream/call_dependent_betabinomial.R）
组织比较通过"组织互为 WT/mutant"的方式复用 dependent 判定框架，跑两次、结果取并集：

| 调用 | --tsv（WT 角色） | genotype1/2（mutant 角色） |
|---|---|---|
| root_vs_flower | `WT_Root_sites.Fisher.tsv` 中 Pass==1 且 `^tRNA-` 的行 | Flower |
| leaf_vs_flower | `WT_Leaf_sites.Fisher.tsv` 中 Pass==1 且 `^tRNA-` 的行 | Flower |

- mutant（Flower）侧文件 = `tRNA/STAT/adjustGap/ncRNA/{Treat,Mock}-WT-Flower.rep{1,2,3}.adjusted.tsv`。
  只用 rep1–3 与 Fisher 表的三重复口径对齐（rep4 是 count-only 补做批次，不引入组织比较）。
- `--genotype2` 与 genotype1 传同一组文件：脚本接口要求两个株系，此处第二次调用即占位，
  AND 判定退化为单比较本身。
- 检验方向：单侧 loss P 已内嵌方向 —— 只有 Flower（mutant 角色）修饰相对参照组织下降
  才可能拿到小 P；delta = 参照组织 Ratio − Flower Ratio。
- 门槛：`--significance q`（BH FDR 家族 = 该次调用的全部 tRNA 候选位点）、`--site-test lrt`、
  `--alpha 0.05`、`--min-delta 0.05`（与其他 dependent 位点分析的效应量门槛一致；
  老板 2026-11-03 拍板，不再沿袭 archive tissue_dependent 的组织比较阈值 0.10）、
  `--min-depth 10`（参照组织 WT 每重复 + Flower Mock/Treat 每重复均须达标）。

**Flower-low 判定 = root_vs_flower OR leaf_vs_flower 任一比较 is_dependent**（并集口径，
与 archive `tissue_dependent_flower_low` 的 pairwise OR 一致）。

## 展示
- 行 = 核 tRNA 家族（family_top 中 `tRNA-Nu-` 成员，代表成员为该家族 mean_CPM 最高者；
  行序：氨基酸字母序 → 反密码子字母序，eMet/iMet 按 Met 归组）。
- 列 = Sprinzl canonical 0–76 位，按 17 个结构区段分块。
- 格子颜色 = 该家族该 canonical 位上 Flower-low 位点的 Flower 修饰比例
  （`WT_Flower_sites.Fisher.tsv` 的 `Treat_avg_ratio`，多拷贝/多位点按 `Treat_avg_depth` 加权合并）；
  白格 = 无 Flower-low 位点。非核 tRNA（Mt/Chl）与 Sprinzl 坐标缺失的位点不入图，单独落表。
- Sprinzl 坐标取自三组织 `WT_*_sites.positive.Sprinzl.tsv` 的并集（同一位点的注释与组织无关）。

## 输入
- `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.Fisher.tsv`（最新一轮分析）
- `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv`
- `tRNA/STAT/adjustGap/ncRNA/{Treat,Mock}-WT-Flower.rep{1,2,3}.adjusted.tsv`
- `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`
- `BID-transcript.script/downstream/call_dependent_betabinomial.R`

## 输出
- `betabinomial/root_vs_flower/`、`betabinomial/leaf_vs_flower/` — bb 脚本原始输出
  （line-level 全位点表 + dependent.tsv + summary.json）
- `flower_low_sites.tsv` — OR 合并后的 Flower-low tRNA 位点表（两次比较的
  ratio/delta/p/q、Flower 深度、Sprinzl_pos、comparison 来源）
- `flower_low_mapped_sites.tsv` — 入图的核 tRNA 位点（家族折叠前）
- `flower_low_excluded_sites.tsv` — 不入图的位点（非核 tRNA / family_top 外 / Sprinzl 缺失）
- `flower_low_trna_heatmap.{pdf,png,svg}` — 展示图

## 运行
```bash
./data.sh            # 两次 bb 调用（并行）→ merge_flower_low.R → plot.R
```
