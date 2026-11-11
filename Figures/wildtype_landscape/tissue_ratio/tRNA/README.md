# tissue_ratio/tRNA — tRNA sites with significantly reduced modification fraction in Flower

## Purpose
Compare modification fractions of tRNAs (`^tRNA-` features in the ncRNA reference) among wild-type Root, Leaf, and Flower tissues using a beta-binomial statistical framework to identify sites with **significantly reduced modification fraction in Flower**, projected onto representative nuclear tRNA families across canonical positions 0–76.

## Methods (BID-transcript.script/downstream/call_dependent_betabinomial.R)
Tissue comparisons repurpose the dependent-calling framework by setting pairwise tissue comparisons as reciprocal reference/mutant contrasts, executed in parallel and combined via union:

| Comparison | --tsv (Reference role) | genotype1/2 (Mutant role) |
|---|---|---|
| root_vs_flower | Rows with Pass==1 and `^tRNA-` in `WT_Root_sites.Fisher.tsv` | Flower |
| leaf_vs_flower | Rows with Pass==1 and `^tRNA-` in `WT_Leaf_sites.Fisher.tsv` | Flower |

- Mutant (Flower) libraries: `tRNA/STAT/adjustGap/ncRNA/{Treat,Mock}-WT-Flower.rep{1,2,3}.adjusted.tsv`. Replicates 1–3 are used to match the triplicates in the Fisher tables (rep4 is a count-only supplementary batch excluded from tissue comparisons).
- `--genotype2` receives the identical set of files as `--genotype1`: the script interface expects two strains, so this duplicate input acts as a placeholder, causing the AND filter to evaluate as a single pairwise test.
- Test direction: The one-sided loss test intrinsically tests for reduction; only sites with decreased modification in Flower relative to the reference tissue yield small p-values; delta = Reference Ratio − Flower Ratio.
- Cutoffs: `--significance q` (BH FDR applied across all tRNA candidates in each run), `--site-test lrt`, `--alpha 0.05`, `--min-delta 0.05` (consistent with effect size thresholds across dependent site analyses), `--min-depth 10` (required across every replicate of reference WT and Flower Mock/Treat).

**Flower-low determination: is_dependent in either root_vs_flower OR leaf_vs_flower** (union approach).

## Visualization
- Rows: Nuclear tRNA families (`tRNA-Nu-` members in `family_top`; representative member is the highest `mean_CPM` transcript; ordered by amino acid alphabetical order → anticodon alphabetical order, with eMet/iMet grouped under Met).
- Columns: Sprinzl canonical positions 0–76, partitioned into 17 structural segments.
- Cell color: Flower modification fraction (`Treat_avg_ratio` from `WT_Flower_sites.Fisher.tsv`, merged across multiple gene copies or sites using `Treat_avg_depth`-weighted averaging) for Flower-low sites. White cells indicate positions without Flower-low sites. Non-nuclear tRNAs (Mt/Chl) and sites lacking Sprinzl coordinates are excluded from the heatmap and output to separate tables.
- Sprinzl coordinates are pooled from the union of `WT_*_sites.positive.Sprinzl.tsv` across the three tissues (annotations per site are tissue-independent).

## Inputs
- `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.Fisher.tsv`
- `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv`
- `tRNA/STAT/adjustGap/ncRNA/{Treat,Mock}-WT-Flower.rep{1,2,3}.adjusted.tsv`
- `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`
- `BID-transcript.script/downstream/call_dependent_betabinomial.R`

## Outputs
- `betabinomial/root_vs_flower/`, `betabinomial/leaf_vs_flower/`: Raw beta-binomial outputs (line-level candidate tables, `dependent.tsv`, `summary.json`)
- `flower_low_sites.tsv`: Union table of Flower-low tRNA sites (ratio, delta, p, q, Flower depth, Sprinzl_pos, source comparison)
- `flower_low_mapped_sites.tsv`: Nuclear tRNA sites included in the heatmap (pre-family collapse)
- `flower_low_excluded_sites.tsv`: Excluded sites (non-nuclear tRNA, non-family_top, or missing Sprinzl coordinates)
- `flower_low_trna_heatmap.{pdf,png,svg}`: Heatmap figure

## Run
```bash
./data.sh            # Two parallel beta-binomial calls → merge_flower_low.R → plot.R
```
