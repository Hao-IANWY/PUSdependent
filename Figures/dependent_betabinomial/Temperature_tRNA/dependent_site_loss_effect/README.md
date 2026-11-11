# dependent_site_loss_effect — Scatter plots of modification loss at TRUA1B-dependent Ψ sites (tRNA)

## Purpose

These scatter plots illustrate the reduction in tRNA Ψ site modification fraction caused by TRUA1B deficiency at 21 °C and 27 °C. X-axis: mutant modification fraction; Y-axis: WT modification fraction. Grey dots represent all WT-positive sites at the given temperature; red dots represent identified TRUA1B-dependent sites. Red dots below the diagonal (WT > Mut) indicate TRUA1B-dependent modification loss.

The mRNA / rRNA version is located in the sibling directory `Figures/dependent_betabinomial/Temperature_mRNA/dependent_site_loss_effect/`.

## Directory layout

```
Figures/dependent_betabinomial/Temperature_tRNA/dependent_site_loss_effect/
├── data_prep.sh                              Logic for generating the two input TSV files
├── plot.R                                    Script to generate PDF/PNG/SVG figures
└── tRNA_<temp>_dependent_site_loss_effect.tsv / .pdf / .png / .svg
```

## Data sources

Upstream directory: `Temperature_tRNA/STAT/dependent_betabinomial/ncRNA/<temp>/trua1b/`. Genotypes include `trua1b_1` and `trua1b_2`. Each genotype × temperature × treatment group has three biological replicates (min_pass_reps = 3). Inputs are based on the beta-binomial LRT calling criteria (as defined in `Temperature_tRNA/STAT/dependent_betabinomial/run.sh`).

Each combination uses three tables:
- `trua1b_1_<temp>.dependent-sites.tsv` — Site-wise table for genotype 1, including all WT candidate sites.
- `trua1b_2_<temp>.dependent-sites.tsv` — Site-wise table for genotype 2, including all WT candidate sites.
- `trua1b-<temp>.dependent.tsv` — Final TRUA1B-dependent sites (intersection of both genotypes).

### Calling criteria

Beta-binomial LRT; BH q < 0.05 (positive direction, meth.diff > 0); effect size `Delta_Ratio(WT-Mut) >= 0.05`; Mock mean depth >= 10. Only sites meeting these criteria in both genotypes (intersection) are classified as dependent.

### tRNA ^tRNA- filtering

The reference `tRNA-unique_sn_snoRNA.fa` contains snRNA and snoRNA sequences, which must be filtered out using the `^tRNA-` prefix (consistent with `Figures/dependent_betabinomial/tRNA` and the sibling `temperature_responsive_scatter` directory). Upstream tables hold 547 (21 °C) / 555 (27 °C) candidates; 493 remain after filtering at both temperatures.

### Site counts

| RNA type | 21 °C WT sites | 21 °C dependent | 27 °C WT sites | 27 °C dependent |
|---|---|---|---|---|
| tRNA | 493 | 0 | 493 | 0 |

No TRUA1B-dependent sites are detected at either temperature. Plots are generated as placeholders to maintain consistency; `plot.R` logs a message when a figure contains zero dependent sites to distinguish an intentional empty panel from a missing data path.

## Usage

```bash
bash data_prep.sh     # Generate 2 TSV files
Rscript plot.R        # Generate 2 figures
```

## Outputs

Two sets of `tRNA_<temp>_dependent_site_loss_effect.{pdf,png,svg}` files. Dimensions: 6.8 × 6.8 in, 600 dpi. Includes `theme_base` border, legend in the bottom-right corner, and temperature annotation in the top-left corner.
