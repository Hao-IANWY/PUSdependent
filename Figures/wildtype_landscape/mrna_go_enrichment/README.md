# mrna_go_enrichment

## 目的
基于当前野生型三个组织的 mRNA 阳性位点结果（Root、Leaf、Flower），按位点的组织共有/特异性进行 GO 富集分析，并使用 `sci-plot` 的 GO enrichment 模版绘制 offset-axis 横向条形图。

组合图 `mrna_go_enrichment_bar_combined` 只画 Roots / Leaves / Flowers 三组，每组为该组织**所有**带修饰 mRNA（基因，即该组织全部阳性位点所属基因的并集）的富集结果。

## 分类规则
- `root_unique`：同一位点只出现在 Root 阳性结果中。
- `leaf_unique`：同一位点只出现在 Leaf 阳性结果中。
- `flower_unique`：同一位点只出现在 Flower 阳性结果中。
- `common`：同一位点在三个组织中均出现。
- 仅在两个组织中出现的位点不进入上述任一 GO 类别，只在 `mrna_go_input_summary.tsv` 中作为 pairwise-only QC 计数记录。
- 转录本 ID 去除 `.` 后版本号后作为 TAIR gene ID，用于 `org.At.tair.db` / `clusterProfiler::enrichGO`。

## sci-plot 运行流程
```bash
bash data_prep.sh
Rscript plot.R
```

`data.sh` 是兼容入口，仅调用 `data_prep.sh`：

```bash
bash data.sh
```

如需复现论文归档中的旧版 dot plot，可运行：

```bash
Rscript plot_mrna_go_enrichment.R
```

## 主要输入
- `../../../polyARNA/STAT/fisher/mRNA/WT-Root_sites.positive.tsv`
- `../../../polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv`
- `../../../polyARNA/STAT/fisher/mRNA/WT-Flower_sites.positive.tsv`

## 数据接口
`data_prep.sh` 生成七个单列、无 header 的 TAIR gene-ID TSV，作为 `plot.R` 的唯一输入接口：

- `root_unique_gene_list.tsv`
- `leaf_unique_gene_list.tsv`
- `flower_unique_gene_list.tsv`
- `common_gene_list.tsv`
- `root_all_gene_list.tsv`（Root 全部带修饰基因，组合图 Roots 组）
- `leaf_all_gene_list.tsv`（Leaf 全部带修饰基因，组合图 Leaves 组）
- `flower_all_gene_list.tsv`（Flower 全部带修饰基因，组合图 Flowers 组）

当前输入基因数：
- Root-unique：154
- Leaf-unique：145
- Flower-unique：60
- Common：104
- Roots（全部）：284
- Leaves（全部）：295
- Flowers（全部）：202

其他中间文件：
- `root_unique_sites.tsv`
- `leaf_unique_sites.tsv`
- `flower_unique_sites.tsv`
- `common_sites.tsv`
- `mrna_go_input_summary.tsv`

## GO 富集设置
- Ontology：BP / CC / MF
- ID conversion：TAIR -> ENTREZID
- Multiple-testing correction：BH
- Significance cutoffs：raw p < 0.05 且 q < 0.05
- Term redundancy reduction：`clusterProfiler::simplify(cutoff = 0.7, by = "p.adjust")`
- 图中每个组织、每个 ontology 展示显著性最强的前 5 个条目（不足 5 个时全部展示）。

## 模版式输出
每组均输出 PDF / PNG / SVG 三种格式：

- `Root_unique_GO_enrichment_bar.pdf/png/svg`
- `Leaf_unique_GO_enrichment_bar.pdf/png/svg`
- `Flower_unique_GO_enrichment_bar.pdf/png/svg`
- `Common_GO_enrichment_bar.pdf/png/svg`
- `Roots_GO_enrichment_bar.pdf/png/svg`
- `Leaves_GO_enrichment_bar.pdf/png/svg`
- `Flowers_GO_enrichment_bar.pdf/png/svg`

组合图（仅 Roots / Leaves / Flowers 三列，每列为该组织全部带修饰 mRNA 的富集结果）：

- `mrna_go_enrichment_bar_combined.pdf/png/svg`

图中：
- x 轴为 `-log10(p.adjust)`；
- 行按 BP / CC / MF 分面；
- 条形颜色区分 ontology；
- GO term 标签绘制在条形起点处；
- 使用 `theme_sciplot_offset()`、`guide_prism_offset()` 和 Arial Bold；
- PDF/SVG 保留可编辑文本。

当前结果中，Flowers（全部 202 个输入基因、182 个映射基因）在 raw p < 0.05 且 BH q < 0.05 下没有显著 GO 条目；脚本会输出空结果表和 “no significant GO terms” 占位图，并在组合图题注中标明。

## 富集结果表
完整结果（clusterProfiler 全部显著条目）：
- `Root_unique_GO_enrichment_bar.csv`
- `Leaf_unique_GO_enrichment_bar.csv`
- `Flower_unique_GO_enrichment_bar.csv`
- `Common_GO_enrichment_bar.csv`
- `Roots_GO_enrichment_bar.csv`
- `Leaves_GO_enrichment_bar.csv`
- `Flowers_GO_enrichment_bar.csv`
- `mrna_go_enrichment_bar_results.csv`（与组合图对应，仅含 Roots / Leaves / Flowers 三组）

## 保留的旧版结果
以下文件来自旧版 dot-plot 工作流，仅作对照保留：
- `*_GO.tsv`
- `*_go_enrichment.pdf/png/svg`
- `mrna_go_enrichment_combined.pdf/png/svg`
- `mrna_go_enrichment_summary.tsv`
