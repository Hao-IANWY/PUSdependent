# representative_nuclear_trna_heatmap — Nuclear tRNA positive site modification ratios across canonical secondary structures

## Purpose
Panel E ("tRNA landscape"): Displays positive-site modification ratios of representative nuclear tRNAs across the **canonical tRNA secondary structure layout** for each of the three wild-type tissues (Root as primary panel; Leaf and Flower generated simultaneously by the same script).

- **Rows** = 46 nuclear tRNA families (representative member = tRNA with the highest `mean_CPM` in `family_top`).
- **Columns** = Canonical positions 0–76 (partitioned into p1, Accstem, Dstem, Dloop, Anticodon stem-loop, Vregion, Tstem, Tloop, CCA; structural color bar, position numbers, and segment names at the bottom).
- **Values** = Depth-weighted mean modification ratio across **all gene copies** for each family at each structural position. Cells without positive sites are shaded light gray (`#FCFCFC`).
- **Color palette** = `#FFF6EB → #C44E00` (Ratio 0–1).

## Representative tRNA and Sprinzl coordinate sources

- **Representative tRNA**: `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv`. Each anticodon family (compartment + amino acid + anticodon) occupies one row, where `feature_id` denotes the member with the highest `mean_CPM` (mean across all Mock libraries). The heatmap displays the **46 nuclear (Nu) families**, ordered alphabetically by `family`; representative members are recorded in the `Representative_ID` column of the mapping table.
- **Sprinzl coordinates**: Taken directly from the `Sprinzl_pos` column in `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv` (annotated by `BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py` using complete annotations from `plantRNA_database_unique.tsv`). This directory performs no coordinate conversions; `plot.R` looks up canonical positions solely to map them onto display columns.
- Each tRNA copy is annotated with its own segment structure via `plantRNA_database_unique.tsv`, eliminating length-mismatch exclusions relative to family representatives.
- Restricted to nuclear tRNAs; mitochondrial (Mt) and plastid (Pt) tRNAs as well as sn/snoRNAs are excluded (recorded as `excluded_nonNu` in the summary table).

## Run

```bash
bash data.sh          # Collapsing (awk) + plotting (Rscript plot.R)
```

Plotting and styling are managed by `Figures/shared/plot_style.R` and `plot_tokens.R` (`save_plot()` outputs PDF, PNG, and SVG simultaneously, standardized to Arial font).

## Input data
- `../../../tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv` (Columns: `chr, pos, Sprinzl_pos, strand, ratio, depth`)
- `../../../tRNA/count/ncRNA/family_top/family_top_tRNA.tsv` (Family representatives defined by expression level)

## Data interface

`data.sh` produces three mapping tables serving as the sole inputs for `plot.R`:

`representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv`, with columns:
`tissue, Rep_ID, Representative_ID, chr, pos, strand, ratio, depth, canonical_pos`

- `Rep_ID`: Family identifier (e.g. `Nu-Arg-TCG`); all copy sites belonging to the same family are collapsed into this row.
- `Representative_ID`: The most highly expressed copy in `family_top` for that family (e.g. `tRNA-Nu-Arg-TCG-2`).
- `canonical_pos`: Annotated `Sprinzl_pos` (0–76).
- When multiple raw sites map to the same (family, canonical position), their ratios are merged using depth-weighted averaging.

Current dataset summary (`trna_landscape_summary.tsv`):

| tissue | sites_in | excluded_nonNu | unmapped | sites_mapped | reps_covered | cells |
|---|---|---|---|---|---|---|
| Root | 439 | 84 | 0 | 355 | 46 | 126 |
| Leaf | 567 | 167 | 0 | 400 | 46 | 134 |
| Flower | 550 | 154 | 0 | 396 | 46 | 132 |

## Outputs
- `representative_nuclear_trna_heatmap_{root,leaf,flower}.pdf/png/svg` (18 in wide × 11.08 in high, 600 dpi)
- `representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv`
- `trna_landscape_summary.tsv`

## Figure-reading notes
- Biological sanity: High-ratio vertical bands appear across the T loop (canonical positions 54–60).
- Typography: PDF embeds Arial-BoldMT subset; SVG specifies `font-family: "Arial"`.

## Known benign warnings
When instantiating heatmap objects without an active graphic device, `ComplexHeatmap::Heatmap()` measures text widths using R's default temporary pdf device. If Arial is not registered in the ghostscript/postscript font database, warnings stating that Arial cannot be found in the PostScript font database may be emitted. Because these metrics match those under cairo devices and a single object is reused across PDF, PNG, and SVG outputs to maintain layout fidelity, `plot.R` selectively silences this specific warning while allowing other warnings to pass through.
