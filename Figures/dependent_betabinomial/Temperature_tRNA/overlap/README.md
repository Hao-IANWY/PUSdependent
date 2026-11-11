# overlap — Venn diagram of 27 °C up-regulated sites vs. 27 °C TRUA1B-dependent sites (tRNA)

## Purpose

This Venn diagram evaluates whether tRNA sites with increased modification at 27 °C overlap with sites requiring TRUA1B for modification at 27 °C. Left circle: 27 °C up-regulated sites (WT_27 vs WT_21, two-sided test, `up_in_WT_27`); Right circle: 27 °C TRUA1B-dependent sites.

The mRNA / rRNA version is located in the sibling directory `Figures/dependent_betabinomial/Temperature_mRNA/overlap/`.

## Data sources

| Set | Upstream directory |
|---|---|
| 27 °C up-regulated | `Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.response.tsv` |
| 27 °C dependent | `Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/27/trua1b/trua1b-27.dependent.tsv` |

Site ID: `chr + pos + strand`. The two groups are combined and marked as `up_only`, `dependent_only`, or `both`. Both sets are filtered by `^tRNA-` to exclude snRNA/snoRNA. Inputs are based on the beta-binomial LRT calling criteria (as defined in `Temperature_tRNA/STAT` processes). Each genotype × temperature × treatment group has three biological replicates (min_pass_reps = 3).

## Results: Overlap is zero

| RNA type | 27 °C up-regulated | 27 °C dependent | Intersection |
|---|---|---|---|
| tRNA | 4 | 0 | 0 |

Four 27 °C up-regulated sites are identified (`tRNA-Nu-Arg-CCG-1:68`, `tRNA-Nu-Leu-TAA-1:65`, `tRNA-Nu-Leu-TAA-2:65`, `tRNA-Pt-His-GTG-1:30`; all `up_only`), but no TRUA1B-dependent sites are detected in tRNA at 27 °C; thus, the intersection is zero — modification gain at 27 °C and TRUA1B dependence define disjoint site sets. Plots are generated as placeholders to maintain consistency; `plot.R` logs a message when a figure contains zero dependent sites.

## Why hand-drawn Venn circles

`plot.R` uses `geom_polygon` to render custom circles. This approach ensures consistent layout and allows empty sets to be represented, whereas standard `VennDiagram` packages fail when intersection counts are zero.

## Usage

```bash
bash data_prep.sh     # Generate overlap TSV
Rscript plot.R        # Generate Venn diagram
```

## Outputs

- `tRNA_27up_vs_trua1b_dependent_overlap.tsv` — Site-wise labeling table
- `tRNA_27up_vs_trua1b_dependent_venn.{pdf,png,svg}` — 9 × 7 in, 600 dpi
