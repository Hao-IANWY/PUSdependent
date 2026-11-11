# merge_heatmap — Cross-genotype union heatmap of dependent tRNA Ψ sites

## Purpose

Per tissue, integrate the dependent tRNA Ψ sites of all 7 genotypes (trua1a/trua1b/trua2/trua3/trua4/trua5/trua6) and display the modification fractions of the union sites in WT and each mutant. Plot styling follows `Figures/shared`.

## Plot design

- **Y axis**: tRNA family (compartment + amino acid + anticodon). Representative members are taken from
  `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv` (highest mean_CPM member per family), grouped by
  compartment (Cy/Mt/Pt). Sites from **all copies of a family** are folded into the same row; non-
  representative copies are never discarded.
- **X axis**: canonical position `fix_pos`; under each position, modification fractions are shown for
  WT, trua1a, trua1b, trua2, trua3, trua4, trua5 and trua6 in that order. Position labels read
  `tRNA_<fix_pos>`.
- **Site set**: union of `../<genotype>/<genotype>_all_tissues.dependent.tsv` over all genotypes within
  the tissue; cells of dependent sites are marked `*`.
- **Modification fractions**: the dependent tables only define the union site set. Fractions for both WT
  and mutants are filled from the all-site ratio tables `Treat-*.adjusted.tsv` under
  `tRNA/STAT/adjustGap/ncRNA/` (aggregated by `data_prep.sh`), so a mutant's fraction is displayed even
  when it differs little from WT. The per-replicate mean of `Treat-WT-<tissue>.rep{1,2,3}.adjusted.tsv`
  matches the Fisher table's `Treat_avg_ratio` site by site; the all-site adjusted tables are used
  because the Fisher tables are filtered subsets.
- **Color scale**: modification fraction 0–1; missing values light gray (`#FCFCFC`).

## Replicate handling

Replicates follow the site-calling exactly: WT rep4 is a count-only replicate added across sequencing
batches and is excluded from calling and from this figure. Mutant fractions are the mean over each
allele's three replicates, then averaged over the two alleles — the same computation chain as the
`<g>_1`/`<g>_2`/`Mut_Ratio` columns of the dependent tables.

## Inputs

- `../trua*/trua*_all_tissues.dependent.tsv` (run each genotype directory's `data_prep.sh` first)
- `tRNA/STAT/adjustGap/ncRNA/Treat-{WT,<genotype>_<line>}-<tissue>.rep{1..3}.adjusted.tsv`
- `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv` (family representatives, expression-based)

## Outputs

- `root/leaf/flower_merged_dependent_site_ratio_heatmap.{pdf,png,svg}`
- matching `_long.tsv` (site × sample long table) and `_matrix.tsv` (fraction matrix)
- `merged_sample_ratios.tsv` (the `data_prep.sh` aggregate; sole fraction input of `plot.R`)

## Usage

```bash
bash data_prep.sh   # aggregate adjusted.tsv → merged_sample_ratios.tsv (~1.5M rows)
Rscript plot.R      # render the three heatmaps
```

## Data summary

Union sites per tissue: Root 517 (47 families × 11 positions), leaf 819 (63 × 13), Flower 696
(58 × 12); folding by family discards no copies.
