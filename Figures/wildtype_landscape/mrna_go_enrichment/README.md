# mrna_go_enrichment — Tissue-specific and shared mRNA GO enrichment analysis

## Purpose
Perform Gene Ontology (GO) enrichment analysis on mRNA positive sites across wild-type tissues (Root, Leaf, Flower) partitioned by tissue sharing and specificity, and render horizontal offset-axis bar plots using the `sci-plot` GO enrichment template.

The combined figure (`mrna_go_enrichment_bar_combined`) presents Roots, Leaves, and Flowers panels, where each group represents the enrichment of **all** modified mRNAs (genes; the union of genes harboring positive sites in that tissue).

## Classification rules
- `root_unique`: Sites detected exclusively in Root positive results.
- `leaf_unique`: Sites detected exclusively in Leaf positive results.
- `flower_unique`: Sites detected exclusively in Flower positive results.
- `common`: Sites detected in all three tissues.
- Sites occurring in exactly two tissues do not enter any tissue-specific or common GO category; they are recorded solely as pairwise-only QC counts in `mrna_go_input_summary.tsv`.
- Transcript IDs are stripped of version suffixes (after `.`) to obtain TAIR gene IDs for `org.At.tair.db` and `clusterProfiler::enrichGO`.

## Run
```bash
bash data_prep.sh
Rscript plot.R
```

`data.sh` provides a compatibility entrypoint calling `data_prep.sh`:

```bash
bash data.sh
```

To reproduce earlier dot plots:

```bash
Rscript plot_mrna_go_enrichment.R
```

## Data sources
- `../../../polyARNA/STAT/fisher/mRNA/WT-Root_sites.positive.tsv`
- `../../../polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv`
- `../../../polyARNA/STAT/fisher/mRNA/WT-Flower_sites.positive.tsv`

## Data interface
`data_prep.sh` generates seven single-column, headerless TAIR gene-ID TSV files serving as the exclusive input interface for `plot.R`:

- `root_unique_gene_list.tsv`
- `leaf_unique_gene_list.tsv`
- `flower_unique_gene_list.tsv`
- `common_gene_list.tsv`
- `root_all_gene_list.tsv` (All modified genes in Root, corresponding to the Roots panel in the combined figure)
- `leaf_all_gene_list.tsv` (All modified genes in Leaf, corresponding to the Leaves panel in the combined figure)
- `flower_all_gene_list.tsv` (All modified genes in Flower, corresponding to the Flowers panel in the combined figure)

Current input gene counts:
- Root-unique: 154
- Leaf-unique: 145
- Flower-unique: 60
- Common: 104
- Roots (all): 284
- Leaves (all): 295
- Flowers (all): 202

Intermediate files:
- `root_unique_sites.tsv`
- `leaf_unique_sites.tsv`
- `flower_unique_sites.tsv`
- `common_sites.tsv`
- `mrna_go_input_summary.tsv`

## GO settings
- Ontologies: BP, CC, MF
- ID conversion: TAIR -> ENTREZID
- Multiple-testing correction: Benjamini-Hochberg (BH)
- Significance cutoffs: raw p < 0.05 and q < 0.05
- Term redundancy reduction: `clusterProfiler::simplify(cutoff = 0.7, by = "p.adjust")`
- Visualization: Top 5 most significant terms per tissue and per ontology (all terms displayed if fewer than 5 meet thresholds).

## Outputs
Each group is rendered in PDF, PNG, and SVG formats:

- `Root_unique_GO_enrichment_bar.pdf/png/svg`
- `Leaf_unique_GO_enrichment_bar.pdf/png/svg`
- `Flower_unique_GO_enrichment_bar.pdf/png/svg`
- `Common_GO_enrichment_bar.pdf/png/svg`
- `Roots_GO_enrichment_bar.pdf/png/svg`
- `Leaves_GO_enrichment_bar.pdf/png/svg`
- `Flowers_GO_enrichment_bar.pdf/png/svg`

Combined figure (Roots / Leaves / Flowers columns, showing all modified mRNAs per tissue):

- `mrna_go_enrichment_bar_combined.pdf/png/svg`

Figure specifications:
- x-axis: `-log10(p.adjust)`
- Row facets: BP, CC, MF
- Bar fill color: Distinguished by ontology
- Term labels: Rendered at the origin of each bar
- Typography and theme: `theme_sciplot_offset()`, `guide_prism_offset()`, Arial Bold
- Vector formats: PDF/SVG maintain editable text.

Note: Flowers (all 202 input genes, 182 mapped genes) yield no significant GO terms under raw p < 0.05 and BH q < 0.05. The script outputs an empty results table, renders a "no significant GO terms" placeholder plot, and notes this in the combined figure caption.

## Enrichment results tables
Full results (all significant clusterProfiler terms):
- `Root_unique_GO_enrichment_bar.csv`
- `Leaf_unique_GO_enrichment_bar.csv`
- `Flower_unique_GO_enrichment_bar.csv`
- `Common_GO_enrichment_bar.csv`
- `Roots_GO_enrichment_bar.csv`
- `Leaves_GO_enrichment_bar.csv`
- `Flowers_GO_enrichment_bar.csv`
- `mrna_go_enrichment_bar_results.csv` (Corresponding to the combined figure, containing Roots / Leaves / Flowers)

## Legacy outputs
Earlier dot-plot renderings of the same enrichment results are retained:
- `*_GO.tsv`
- `*_go_enrichment.pdf/png/svg`
- `mrna_go_enrichment_combined.pdf/png/svg`
- `mrna_go_enrichment_summary.tsv`
