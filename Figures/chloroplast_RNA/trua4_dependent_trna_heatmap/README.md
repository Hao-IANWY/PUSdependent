# trua4_dependent_trna_heatmap — Structural heatmap of TRUA4-dependent chloroplast tRNAs in leaf

## Purpose

Profile TRUA4-dependent chloroplast tRNA Ψ sites in leaf tissue, illustrating their distribution across canonical tRNA secondary structure positions (positions 0–76) and the magnitude of modification fraction reduction.

## Inputs

Inputs are derived using beta-binomial calling criteria (`STAT/dependent_betabinomial`):

- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4-leaf.dependent.Sprinzl.tsv` (standard dependent site table requiring both mutant alleles; includes `Sprinzl_pos`)
- `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/trua4_{1,2}-leaf.dependent-sites.tsv` (`Mut_Ratio` per mutant allele, retained as column records)

## Filtering and Mapping Rules

- Retain chloroplast tRNAs with `chr` starting with `tRNA-Pt-` (excluding nuclear/mitochondrial tRNAs and sn/snoRNAs in the same ncRNA reference)
- Canonical Sprinzl positions are taken directly from `Sprinzl_pos` in the input table
- Column layout corresponds to canonical positions 0–76 across 17 structural regions (77 columns in total: p1, Acceptor stem, D stem-loop, Anticodon stem-loop, Variable, T stem-loop, Acceptor stem, CCA), isomorphic to the heatmaps in `Figures/dependent_betabinomial/tRNA`
- Rows correspond to individual dependent chloroplast tRNAs, sorted in ascending order by canonical position and then alphabetically by tRNA identifier, grouping sites at Ψ38, Ψ39, and Ψ40
- Value = `Delta_Ratio(WT-Mut)`; when multiple observations exist for the same (tRNA, canonical position) pair, the maximum value is displayed

## Data Summary

9 dependent sites, all shared across both mutant alleles (`shared`):

| Canonical position | tRNA |
|---|---|
| 38 | Pt-Pro-TGG-1 |
| 39 | Pt-Cys-GCA-1, Pt-His-GTG-1, Pt-Leu-CAA-1, Pt-Leu-TAA-1, Pt-Ser-GCT-1, Pt-Thr-TGT-1, Pt-Tyr-GTA-1 |
| 40 | Pt-Asn-GTT-1 |

TRUA4 targets exclusively canonical positions Ψ38–40 at the 3' end of the anticodon loop and anticodon stem in chloroplast tRNAs. Under beta-binomial calling criteria (AND logic across both mutant alleles), Pt-Val-TAC at Ψ55 does not meet significance (q = 0.63 / 0.87), whereas Pt-Cys-GCA, Pt-Leu-TAA, and Pt-Ser-GCT are identified as dependent, resulting in 9 sites strictly localized to Ψ38–40.

## Outputs

- `trua4_pt_trna_dependent_union.tsv`: Chloroplast dependent site table (`chr, pos, Sprinzl_pos, strand, context, WT_Ratio, Mut_Ratio, Delta_Ratio(WT-Mut), trua4_1, trua4_2, dependency_support`)
- `trua4_pt_trna_dependent_mapped_sites.tsv`: Mapping table annotated with `display_label`, `canonical_pos`, and `matrix_col_id`
- `trua4_pt_trna_dependent_heatmap.pdf/png/svg`

Plot annotation: Rows represent individual dependent chloroplast tRNAs; columns represent canonical positions 0–76; color scale indicates `Delta Ratio` (`#FFF5EE → #8B0000`), with the upper limit calculated using `y_axis_cover()` to enclose the maximum value; cells without detected sites are colored light gray.

## Usage

```bash
bash data.sh
```
