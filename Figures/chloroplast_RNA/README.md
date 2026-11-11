# chloroplast_RNA — Chloroplast Ψ landscape

## Purpose

Expand the investigation from nuclear-encoded transcripts to chloroplast RNAs: profile the number of Ψ sites and modification fraction distributions across chloroplast rRNA, tRNA, and mRNA in WT leaf tissue, provide orthogonal validation using chloroplast-purified RNA sequencing libraries, and characterize TRUA4-dependent Ψ sites on chloroplast rRNA, tRNA, and mRNA along with their structural distribution on tRNAs.

Plot styling follows `Figures/shared/plot_style.R` + `plot_tokens.R` (Arial font, rendering PDF/PNG/SVG via `save_plot()`).

TRUA4-dependent analysis is performed under the beta-binomial calling criteria (`STAT/dependent_betabinomial`), where TRUA4 dependence is required in both mutant alleles (AND logic), consistent with `Figures/dependent_betabinomial`.

## Directory Structure

```
Figures/chloroplast_RNA/
├── chloroplast_site_overview/              WT leaf chloroplast site counts and modification fraction distribution
├── chloroplast_purified_rna_validation/    Chloroplast-purified RNA site counts and overlap Venn diagrams
├── chloroplast_chr_overview/               Circular ChrC genome plot (with mRNA gene annotations)
├── trua4_dependent_site_overview/          WT vs trua4 chloroplast site scatter plots (dependent sites highlighted in red)
├── trua4_dependent_trna_heatmap/           Structural heatmap of chloroplast TRUA4-dependent tRNA sites
└── README.md
```

## Data Sources

| Data | Path |
|---|---|
| WT leaf rRNA | `polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv` |
| WT leaf mRNA | `polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv` |
| WT leaf tRNA | `tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv` |
| Chloroplast-purified RNA | `Chloroplast/STAT/fisher/{ncRNA,mRNA}/Chloroplast_sites.positive.tsv` |
| trua4 mutant alleles modification fraction / dependent | `polyARNA/STAT/dependent_betabinomial/{rRNA,mRNA}/leaf/trua4/`, `tRNA/STAT/dependent_betabinomial/ncRNA/leaf/trua4/` |
| Chloroplast annotation | `Araport11_fixed.gtf` (`ChrC`) |

### Identification Criteria for Chloroplast Sites

- rRNA: `chr == "Ch23S"`
- tRNA: `chr` starting with `tRNA-Pt-` (excluding nuclear/mitochondrial tRNAs and sn/snoRNAs from the same ncRNA reference)
- mRNA: `chr` starting with `ATCG` (chloroplast-encoded genes)

Canonical Sprinzl positions for tRNAs are taken directly from upstream `Sprinzl_pos`.

*Note: `Chloroplast/` is a symbolic link pointing to the sequencing workflow directory and is not tracked in the repository; derived tables used for plotting are included.*

## Data Summary

| Panel | rRNA | tRNA | mRNA |
|---|---|---|---|
| WT leaf sites (site overview) | 4 | 64 | 18 |
| Chloroplast-purified RNA sites | 3 | 59 | 15 |
| Overlap | 3 | 55 | 11 |
| trua4-dependent sites | 2 | 9 | 0 |

All 9 TRUA4-dependent tRNA sites reside at canonical positions Ψ38–40 (38: Pro-TGG; 39: Cys-GCA, His-GTG, Leu-CAA, Leu-TAA, Ser-GCT, Thr-TGT, Tyr-GTA; 40: Asn-GTT).

## Execution Order

```bash
bash Figures/chloroplast_RNA/chloroplast_site_overview/data.sh
bash Figures/chloroplast_RNA/chloroplast_purified_rna_validation/data.sh  # reads output from the previous directory
bash Figures/chloroplast_RNA/chloroplast_chr_overview/data.sh
bash Figures/chloroplast_RNA/trua4_dependent_site_overview/data.sh
bash Figures/chloroplast_RNA/trua4_dependent_trna_heatmap/data.sh
```

`chloroplast_purified_rna_validation` takes `chloroplast_site_overview_density_sites.tsv` from `chloroplast_site_overview` as the WT leaf reference set and must be run after it. The other directories are mutually independent; the two trua4-dependent directories directly read `STAT/dependent_betabinomial`.
