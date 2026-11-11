# temperature_responsive_scatter — Scatter plots of temperature-responsive Ψ sites in WT (21 °C vs. 27 °C) (tRNA)

## Purpose

These scatter plots show WT tRNA Ψ sites at 21 °C and 27 °C on the same axes: X-axis = 21 °C modification fraction; Y-axis = 27 °C modification fraction. The diagonal represents y=x. Sites below the diagonal are down-regulated at 27 °C (blue); sites above are up-regulated (red); others are grey.

The mRNA / rRNA version is located in the sibling directory `Figures/dependent_betabinomial/Temperature_mRNA/temperature_responsive_scatter/`.

## Directory layout

```
Figures/dependent_betabinomial/Temperature_tRNA/temperature_responsive_scatter/
├── data_prep.sh                              Logic for generating the input TSV file
├── plot.R                                    Script to generate PDF/PNG/SVG figures
└── tRNA_temperature_responsive_scatter.tsv / .pdf / .png / .svg
```

## Data sources

Upstream directory: `Temperature_tRNA/STAT/Temperature_response_betabinomial/ncRNA/WT_27_vs_WT_21.all-sites.tsv`. `all-sites.tsv` is the full candidate table: the union of Pass==1 sites from `WT_21_sites.Fisher.tsv` and `WT_27_sites.Fisher.tsv`, with both temperatures' modification fractions and mean Mock depths recomputed from the adjustGap pileup. Sites positive at only one temperature are deliberately retained — they are strong candidates for modification loss at the other temperature. This directory therefore does not recompute from raw fisher/adjustGap outputs; `data_prep.sh` only normalizes column names and filters rows.

### Calling criteria

Beta-binomial GLMM (global theta + site-wise LRT); BH q < 0.05 (two-sided FDR); effect size `|Delta_Ratio(WT_27-WT_21)| >= 0.05`; Mock/Treat replicate depth >= 10. Sites are categorized as `up-regulated` (up_in_WT_27), `down-regulated` (down_in_WT_27), or `non-responsive`. Each genotype × temperature × treatment group has three biological replicates (min_pass_reps = 3).

### tRNA ^tRNA- filtering

The reference `tRNA-unique_sn_snoRNA.fa` contains snRNA and snoRNA sequences, which must be filtered out using the `^tRNA-` prefix.

### Site counts

| RNA type | Candidate sites | Responsive sites | 27 °C up-regulated | 27 °C down-regulated |
|---|---|---|---|---|
| tRNA | 502 | 4 | 4 | 0 |

All four responsive sites are tRNAs (`tRNA-Nu-Arg-CCG-1:68`, `tRNA-Nu-Leu-TAA-1:65`, `tRNA-Nu-Leu-TAA-2:65`, `tRNA-Pt-His-GTG-1:30`, all up-regulated; smallest q = 8.8e-05). The upstream ncRNA table holds 567 candidates; `^tRNA-` filtering retains 502 (65 sn/snoRNA candidates removed), so the scatter contains 4 red points.

## Usage

```bash
bash data_prep.sh     # Generate TSV file
Rscript plot.R        # Generate figures
```

## Outputs

`tRNA_temperature_responsive_scatter.{pdf,png,svg}`. Dimensions: 6.8 × 6.8 in, 600 dpi. Includes `theme_base` border and legend in the bottom-right corner.
