# mRNA_motif — WT mRNA Fisher 阳性位点窗口 seqlogo（Y 轴 bits）

wildtype_landscape 的 mRNA 阳性位点 ±2 / ±5 / ±10 nt 窗口序列 logo，与模型版
（PUScaller `data/predict_result/polyARNA/PileupNet_RNAFM_v1/mRNA_sites_feature/motif/`）
脚本同构、口径对齐，可直接并排比较。

## 输入口径

- 分子：`polyARNA/STAT/fisher/mRNA/WT-{Root,leaf,Flower}_sites.positive.tsv`
  （Root 235 / Leaf 199 / Flower 163）
- 序列：Fisher 阳性表无 context 列，三个窗口均按转录本 1-based 拼接坐标从
  Araport11 代表转录本 cDNA（FASTA）切片（sense 方向，无需反向互补）；
  清洗规则：大写、U→T、仅完整窗口、无 N
- 本轮全保留（无剔除）；中心碱基三窗口 100% T（QC 见各窗口
  `mrna_motif_center_base_stats.tsv`）

## 运行

```bash
bash run.sh
```

## 目录

```
mRNA_motif/
├── run.sh                        # 总入口
├── extract_window_sequences.py   # 三窗口序列提取（FASTA 切片）
├── plot_mrna_motif_seqlogo.R     # ggseqlogo 四面板（Roots/Leaves/Flowers/All），bits 标度
└── +-{2,5,10}/                   # 每窗口：序列表 + 中心碱基 QC + seqlogo 三格式
```

## 结果要点（2026-11-04）

- 位点 0 位 100% T（满 2 bits）；
- -1 位 T 富集弱于模型版（Leaf ~0.3 / Flower ~0.28 / Root ~0.15 bits；模型版 Leaf ~1.0）——
  模型学到的上游 T 上下文偏好显著强于规则法；
- +2 位 G 弱富集（~0.15–0.2 bits，三组织一致），模型版该位置信号很弱；
- ±5/±10 窗口显示上下游信号均衰减至背景，主要信息集中在 ±2。
