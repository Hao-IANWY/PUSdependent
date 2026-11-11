# dependent_betabinomial — Writer-dependent Ψ site figures (beta-binomial calling)

Figure scripts and processed tables for PUS writer-dependent pseudouridine (Ψ) sites, identified with the beta-binomial GLMM caller [`BID-transcript.script/downstream/call_dependent_betabinomial.R`](../../../BID-transcript.script/downstream/call_dependent_betabinomial.R). All data are derived from the `STAT/dependent_betabinomial` outputs of the corresponding dataset directories.

## Calling criteria

A candidate Ψ site (positive in WT) is called writer-dependent when:

- **Statistical significance**: BH-corrected q < 0.05 (one-sided test for loss in the mutant).
- **Effect size**: ΔRatio(WT − Mut) ≥ 0.05 for mRNA and tRNA (ncRNA), ≥ 0.10 for rRNA.
- **Read depth**: mean Mock depth ≥ 10 in each replicate group.
- **Consistency**: the criteria must hold in **both** independent mutant alleles (AND).

## Directory layout

```
Figures/dependent_betabinomial/
├── tRNA/    Full tRNA group: sites_count / loss_effect / motif /
│            merge_heatmap / count_missed_dependent / per-genotype heatmaps,
│            distribution boxplots and motif seqlogos (trua1a–trua6) / site_specificity
├── mRNA/    poly(A)+ RNA: sites_count / loss_effect / ratio_distribution /
│            region_distribution / heatmap (fraction of highly modified WT sites
│            that are dependent, across genotypes) / trua1a_dependent (motif +
│            secondary-structure analyses)
└── rRNA/    poly(A)+ RNA: sites_count / loss_effect / tissue_specificity
```

## mRNA/trua1a_dependent analysis order

1. `motif/data_prep.sh` → `motif/plot.R`: 5-nt context sequence logos of dependent sites per tissue (Root 28 / leaf 14 / Flower 25 sequences).
2. `secondary_str/site_info/build_case_sites.py`: union of the three tissues' dependent rows into a physical case-site table (67 rows → 43 sites: 28 single-tissue / 6 two-tissue / 9 three-tissue).
3. `secondary_str/extract_windows.sh` → `secondary_str/data_prep.sh` (local folding with ViennaRNA 2.7.0 at 21 °C via `fold_structure.py`) → `plot.R` + `plot_union_paired_profile.R` + `plot_tissue_paired_profile.R` (paired-probability profiles in −20…+60 nt windows around sites, per tissue and cross-tissue union).

## Data summary

Dependent site counts under the current calling criteria:

| Group | Dependent sites |
|---|---|
| tRNA (`^tRNA-` filtered; 21 genotype × tissue combinations) | 544 |
| poly(A)+ mRNA | 53 |
| poly(A)+ rRNA | 9 |

Every figure directory's `data_prep.sh` can regenerate its tables directly from the corresponding `STAT/dependent_betabinomial` outputs.
