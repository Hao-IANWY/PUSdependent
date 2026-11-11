# trua4_dependent_site_overview — WT vs trua4 scatter plots of chloroplast Ψ sites in leaf

## Purpose

Compare modification fractions between WT and `trua4` leaf tissue across all chloroplast rRNA, tRNA, and mRNA Ψ sites, highlighting TRUA4-dependent Ψ sites in red.

## Inputs

Input tables follow beta-binomial calling criteria (`STAT/dependent_betabinomial`), consistent with `Figures/dependent_betabinomial`:

- `polyARNA/STAT/dependent_betabinomial/rRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv`
- `polyARNA/STAT/dependent_betabinomial/rRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv`
- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`
- `polyARNA/STAT/dependent_betabinomial/mRNA/leaf/trua4/WT-leaf_sites.Treat_avg_ratio.tsv`
- `polyARNA/STAT/dependent_betabinomial/mRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv`

The WT tables are WT positive site tables (`Treat_avg_ratio` version of `WT-*_sites.positive.tsv`, with identical site sets and `ratio` values).

## Filtering Rules

- rRNA: `chr == "Ch23S"`
- tRNA: `chr` starting with `tRNA-Pt-` (excluding nuclear/mitochondrial tRNAs and sn/snoRNAs in the same ncRNA reference)
- mRNA: `chr` starting with `ATCG` (chloroplast-encoded genes)

## TRUA4-Dependent Calling Criteria

Consistent with `Figures/dependent_betabinomial`: **`is_dependent` must be TRUE in both mutant alleles** to be classified as dependent and highlighted in red (`site_status = Dependent sites`). Sites satisfying the criteria in only one allele are labeled as `trua4_1_only` or `trua4_2_only` in `dependency_support` and are not highlighted.
`Mut_Ratio` is calculated as the mean of `Mut_Ratio` across both mutant alleles.

## Data Summary

| RNA | WT sites | Dependent sites |
|---|---|---|
| rRNA | 4 | 2 |
| tRNA | 64 | 9 |
| mRNA | 18 | 0 |

An additional tRNA site meets the criteria only in `trua4_2` (`trua4_2_only`) and is excluded from dependent sites under the two-allele requirement (AND logic).

## Outputs

- `trua4_dependent_site_overview_dotplot.tsv`: Merged site table across all three RNA classes (`rna_type, chr, pos, strand, site_id, WT_Ratio, Mut_Ratio, Delta_Ratio, site_status, dependency_support, trua4_1, trua4_2, trua4_1_dep, trua4_2_dep`)
- `trua4_dependent_site_overview_summary.tsv`: Total site counts and dependent site counts per RNA class
- `trua4_dependent_site_overview_dotplot_{rrna,trna,mrna}.pdf/png/svg`

Plot annotation: x-axis = `trua4` mutant mean modification fraction, y-axis = WT mean modification fraction; dashed line indicates y = x; gray dots = all WT sites, red dots = dependent sites; axis range = 0–100%.

## Usage

```bash
bash data.sh
```
