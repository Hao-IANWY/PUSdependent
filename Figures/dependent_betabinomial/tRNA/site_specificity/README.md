# dependent_betabinomial/tRNA/site_specificity — Site-specificity figures by Sprinzl position

## Purpose

For a given canonical (Sprinzl) position: how many **representative tRNAs** carry WT modification
there, and what happens to their modification fractions in a set of target mutants (two alleles each)?
Every position gets its own folder with two figures — a heatmap and a modification-fraction
distribution plot. The current set covers the Flower panels:

- `pos27/` `pos28/` `pos29/` — genotypes trua1a + trua5 (anticodon stem; TRUA1A acts at 27/29, TRUA5 at 28)
- `pos38/` `pos39/` `pos40/` — genotypes trua3 + trua6 (anticodon loop; TRUA3 acts at 38, TRUA6 at 39/40)
- `pos67/` `pos68/` `pos69/` — genotype trua1a (T loop)

**Only nuclear (Nu) tRNAs are included**; Mt/Pt compartments are excluded, the heatmap carries no
compartment annotation, and distribution titles omit tissue/sample-size labels.

## Directory layout

```
site_specificity/
├── build_site_specificity_table.sh   Shared data construction (parameterized by
│                                     tissue × genotypes × sprinzl_pos × compartments)
├── site_specificity_common.R         Shared heatmap + distribution implementation
│                                     (per genotype: two fraction columns + two per-line
│                                     dependency columns)
├── WT-Flower_sites.Sprinzl.tsv       Shared intermediate: Sprinzl-annotated WT positive
│                                     universe (shared by all position folders of a tissue)
├── pos27/ … pos69/
│   ├── data_prep.sh                  Binds TISSUE=Flower / GENOTYPES / POS=N / COMPARTMENTS=Nu
│   ├── plot.R                        Declares parameters and selects the two figures
│   ├── pos<N>_flower_site_table.tsv  Selected-site wide table (one row per tRNA site)
│   └── pos<N>_flower_site_specificity_{heatmap,distribution}.{pdf,png,svg}
```

Extending: to add a position, copy any folder and replace the position number (adjust `GENOTYPES` if
needed). To add a tissue, copy a folder and change `TISSUE` and the file names (the WT Sprinzl
intermediate is generated per tissue at the shared level). To add a genotype, extend `GENOTYPES` /
`genotypes` in both scripts and add a light→dark allele color pair to `site_specificity_line_palette`
in `site_specificity_common.R` (existing: trua1a pink/red, trua5 light/dark blue, trua3 light/dark
green, trua6 light/dark purple; a missing entry is a hard error).

## Site inclusion criteria

1. **WT modified** = the site belongs to the WT Fisher positive universe
   (`tRNA/STAT/dependent_betabinomial/ncRNA/<tissue>/<g>/WT-<tissue>_sites.Treat_avg_ratio.tsv`,
   byte-identical across genotype directories and verified with `cmp` during construction; every row
   has Pass = 1, i.e. all three WT Treat replicates show the CMC effect). WT fractions are the table's
   `Treat_avg_ratio`.
2. **Position** = `Sprinzl_pos` equals the target position (native positions differ across copies
   because D-loop and variable-region lengths vary, e.g. position 27 is native 26/27/28/29 in different
   copies). The WT universe itself is not Sprinzl-annotated upstream, so `data_prep.sh` calls
   `add_sprinzl_pos.py` once per tissue to build the annotated shared intermediate.
3. **Representative tRNA** = `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv` (highest mean_CPM member
   per family); `^tRNA-` filters sn/snoRNAs.
4. **Compartment** = the `COMPARTMENTS` parameter (comma list of Nu/Mt/Pt; empty = all); currently Nu
   only. Mt/Pt sites stay out of table and figure.

Per-allele fractions are the `Mut_Ratio` of `<g>_<line>-<tissue>.dependent-sites.tsv` (mean of that
allele's Treat replicates over the unfiltered pileup), consistent with the genotype figure group; a
site missing from an allele's table stays **NA and is never treated as 0% modification**. `data_prep.sh`
additionally verifies site by site that the row-table `WT_Ratio` equals the WT universe's
`Treat_avg_ratio` and that both `context` values agree, and aborts otherwise.

## How to read the figures

**Heatmap** `pos<N>_flower_site_specificity_heatmap` (styling aligned with `../../mRNA/heatmap`)

- Rows = nuclear representative tRNAs modified at the position in WT, ordered by WT fraction
  (descending); row labels are family names (`tRNA-Nu-` prefix and copy number stripped, e.g. `Arg-CCT`).
- Columns = WT plus both alleles of each genotype; cell values are modification fractions (one decimal
  below 10%, integer otherwise).
- Red frames mark **per-line** beta-binomial dependency (a single qualifying allele is framed, so the
  framed count can exceed the both-allele count in the table below by 1); gray NA cells mean no
  measurement for that allele.

**Modification-fraction distribution** `pos<N>_flower_site_specificity_distribution`
(styling aligned with the `../trua1a` position-distribution boxplots + jittered points)

- x = the five sample groups, y = Ψ modification fraction (%); single panel titled "Ψ<N> WT-modified tRNAs".
- WT gray; each genotype one light→dark allele pair (trua1a pink, trua5 blue, trua3 green, trua6 purple).

## Data summary (Flower, Nu-only)

| Position | tRNAs | trua1a | trua5 | trua3 | trua6 | Median WT → loss (writer) |
|---|---|---|---|---|---|---|
| Ψ27 | 13 | **11** | 0 | — | — | 46.9% → 5.6/4.9% (1a) |
| Ψ28 | 4 | 1 | **4** | — | — | 40.6% → 12.8/18.4% (5) |
| Ψ29 | 1 (Val-TAC) | **1** | 0 | — | — | 43.4% → 4.6/4.1% (1a) |
| Ψ38 | 5 | — | — | **5** | 0 | 26.9% → 7.3/4.7% (3) |
| Ψ39 | 16 | — | — | 1 | **12** | 61.8% → 13.3/14.9% (6) |
| Ψ40 | 1 (Ser-GCT) | — | — | 0 | **1** | 73.9% → 15.0/18.8% (6) |
| Ψ67 | 1 | **1** | — | — | — | 51.8% → 6.1/3.0% (1a) |
| Ψ68 | 6 | **4** | — | — | — | 41.8% → 5.7/3.6% (1a) |
| Ψ69 | 1 (Ser-CGA) | **1** | — | — | — | 19.7% → 3.2/1.8% (1a) |

(Bold = sites dependent in both alleles simultaneously; "—" = genotype not part of that position's panel.)

The two writer groups show mirror-image specificity: Ψ27/29 belong to TRUA1A and Ψ28 to TRUA5
(anticodon stem); Ψ38 belongs to TRUA3 and Ψ39/40 to TRUA6 (anticodon loop); T-loop Ψ67/68/69 also
belong to TRUA1A. Notable single sites:

- tRNA-Nu-Cys-GCA-10: Ψ27 is unaffected in trua1a, but its Ψ28 is dependent in both trua1a and trua5 —
  two positions on one tRNA fall to different writers with crosstalk.
- The trua1a exception at Ψ27, Arg-TCG-2 (0.23 → 0.06): the fraction drops but the per-line q misses
  the threshold.
- Ψ39 Trp-CCA-1 qualifies only in trua6_1 (trua6_2 FALSE) — the framed-count vs both-TRUE difference.
- Ψ68 Arg-CCG-1 and Leu-AAG-2 qualify in neither allele of trua1a (frames 8 = 4×2).

## Usage

```bash
cd pos27          # or any other position folder
bash data_prep.sh # build tables only, no plotting
Rscript plot.R    # render figures only
```
