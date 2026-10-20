# call_dependent_methylKit.R 统计学检验原理

> 脚本位置：`BID-transcript.script/downstream/call_dependent_methylKit.R`
> 配套示意图：`dependent_methylKit_test_schematic.pdf/.png/.svg`（`plot.R` 生成）

## 一句话概括

对每一个候选位点，把「WT 3 个重复 vs 突变株 3 个重复」的修饰计数摆在一起，
用 **带覆盖度加权、带过离散校正的 logistic 回归** 回答一个问题：
**突变组的修饰率是不是显著比 WT 低？** 最后再叠加效应量门槛和两个株系的交集。

## 第一步：把 ratio 还原成"命中次数"的计数

每个样本在每个位点上的原始数据是 `(depth, ratio)`。脚本把它换算回甲基化式的计数：

- 命中数（被修饰的 U 数）= `round(depth × ratio)`
- 未命中数 = `depth − 命中数`

统计检验**不用 ratio 本身，而用这对整数计数**。这样覆盖度高的位点天然话语权更大：
10/20 和 1/2 虽然都是 50%，但前者证据强得多，回归会自动体现这一点。

## 第二步：零覆盖重复按缺失处理（预处理）

adjusted 文件是全量位点网格，没测到的位点被硬写成 `depth=0, ratio=0`——
这不是"修饰丢了"，是"没数据"。对策：

- 算 Mut_Ratio 均值时，depth=0 的重复按**缺失**跳过；
- 检验前把位点按"哪些重复是零覆盖"分成不同模式，**剔除对应样本后单独检验**
  （例如突变株 rep3 没覆盖，该位点就用 3v2 检验），位点本身保留；
- 若某一组样本被剔光，则无法检验，p 留 NA（status = "failed"）。

这与 `call_dependent_stats.R` 中 0-trial 观测在 GLM 里自动零权重的效果对齐。

## 第三步：核心模型——加权二项 logistic 回归

本项目每组 3 个生物学重复（3v3），methylKit 走 **logReg** 路径。
对每个位点独立拟合：

```
logit(π) = β0 + β1 · x
```

- `π`：该位点的修饰概率；
- `x`：分组标签，**WT = 1，突变株 = 0**；
- 每条观测按 depth（coverage）加权；
- **β1 就是"WT 比突变株高多少"的量化**：β1 > 0 → WT 修饰更高，
  正是"依赖位点（突变株修饰丢失）"的方向。

（只有当总共 2 个样本、无重复时，methylKit 才强制退回 Fisher 精确检验
fast.fisher——本项目 3v3 用不到。）

## 第四步：过离散校正——F 检验

生物学重复间的变异往往比理想二项分布大（overdispersion）。methylKit 用
`overdispersion = "MN"` 估计过离散因子 φ：

- **φ > 1**（重复间确实更散）→ 用 `deviance/φ ~ F(1, n−2)` 做 **F 检验**。
  相当于把"重复不齐"算进误差里，p 值变保守，防止假阳性；
- **φ ≤ 1**（不比二项更散）→ 退化为普通**卡方检验**。

这就是 `--methylkit-test F` 的含义，也是它和普通 Fisher 精确检验的关键区别：
**Fisher 只看合并计数、看不见重复内变异；logistic + F 显式利用了 3 个重复的信息。**

检验输出的是**双侧 p 值**（"两组有差异"）。

## 第五步：BH 多重校正

methylKit 的 `calculateDiffMeth(adjust = "BH")` 会对同一比较的全部候选位点做
Benjamini–Hochberg 校正，把双侧 p 折算成 q 值：

```
q = m · p / rank
```

q 的含义是"在我叫出的这批位点里，假阳性的预期比例"（FDR）。BH 不估计真
null 比例（等价于假设 π0 = 1），因此在候选位点少、真信号多时会比 SLIM 保守。

## 第六步：方向 + 判定——一个位点要"通关"需同时满足 4 个条件

| # | 条件 | 含义 |
|---|------|------|
| 1 | Mock 平均 depth ≥ 10 | 突变株 mock 样本里这个位点测得够深，"丢失"才可信 |
| 2 | BH q < 0.05 且 meth.diff > 0 | WT 显著高于突变株（FDR 校正后，方向必须为正） |
| 3 | Δ = WT_Ratio − Mut_Ratio ≥ 0.05 | 效应量门槛，挡住"显著但只差 1%"的位点 |
| 4 | **两个株系各自都满足 1–3** | AND 合并取交集（如 trua1a_1 ∩ trua1a_2） |

注意口径分工：q 值用 methylKit 的加权模型算（双侧 FDR），方向由 meth.diff 判断；
Δ 用的是**未筛选的简单均值差**（WT_Ratio 与突变株 3 个 Treat 重复 ratio 的简单均值），
meth.diff 不参与 Δ 计算。

## 流程总览

```
候选位点
  → 每样本 (depth, ratio) 还原成计数 (k, depth−k)
  → depth=0 重复按缺失剔除（3v3 可退化为 3v2）
  → 加权 logistic 回归：logit(π) = β0 + β1·x
  → φ 校正：φ>1 走 F 检验，φ≤1 走卡方
  → 双侧 p
  → BH q（同一比较全部候选位点，双侧 FDR）
  → mock depth ≥ 10 且 q < 0.05 且 meth.diff > 0 且 Δ ≥ 0.05
  → 两个株系取交集
  = dependent 位点
```

## 参考

- Akalin et al. 2012, *Genome Biology* 13:R87（methylKit 经典工作流）
- 同项目并列实现：`call_dependent_stats.R`（quasi-binomial / glmmTMB 计数模型版）
