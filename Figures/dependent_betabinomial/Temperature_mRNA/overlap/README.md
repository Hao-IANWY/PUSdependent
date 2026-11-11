# overlap — Venn diagram of 27 °C up-regulated sites vs. 27 °C TRUA1B-dependent sites (mRNA / rRNA)

## Purpose

This Venn diagram evaluates whether sites with increased modification at 27 °C overlap with sites requiring TRUA1B for modification at 27 °C. Left circle: 27 °C up-regulated sites (WT_27 vs WT_21, two-sided test, `up_in_WT_27`); Right circle: 27 °C TRUA1B-dependent sites.

The tRNA version is located in the sibling directory `Figures/dependent_betabinomial/Temperature_tRNA/overlap/`.

## Data sources

| Set | Upstream directory |
|---|---|
| 27 °C up-regulated | `Temperature/STAT/Temperature_response_betabinomial/<type>/WT_27_vs_WT_21.response.tsv` |
| 27 °C dependent | `Temperature/STAT/dependent_betabinomial/<type>/27/trua1b/trua1b-27.dependent.tsv` |

Site ID: `chr + pos + strand`. The two groups are combined and marked as `up_only`, `dependent_only`, or `both`. Inputs are based on the beta-binomial LRT calling criteria (as defined in `Temperature/STAT/Temperature_response_betabinomial` and `dependent_betabinomial`). Each genotype × temperature × treatment group has three biological replicates (min_pass_reps = 3).

## Results: No overlap

| RNA type | 27 °C up-regulated | 27 °C dependent | Intersection |
|---|---|---|---|
| mRNA | 0 | 0 | 0 |
| rRNA | 0 | 0 | 0 |

No sites were identified as 27 °C up-regulated or TRUA1B-dependent in either RNA class; modification gain at 27 °C and TRUA1B dependence define disjoint site sets. Plots are generated as placeholders to maintain consistency; `plot.R` logs a message when a figure contains zero dependent sites.

## Why hand-drawn Venn circles

`plot.R` uses `geom_polygon` to render custom circles. This approach ensures consistent layout and allows empty sets to be represented, whereas standard `VennDiagram` packages fail when intersection counts are zero.

## Usage

```bash
bash data_prep.sh     # Generate 2 overlap TSV files
Rscript plot.R        # Generate 2 Venn diagrams
```

## Outputs

- `<rna>_27up_vs_trua1b_dependent_overlap.tsv` — Site-wise labeling table
- `<rna>_27up_vs_trua1b_dependent_venn.{pdf,png,svg}` — 9 × 7 in, 600 dpi
