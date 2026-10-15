# tissue_ratio

## Purpose
Compares WT pseudouridine modification ratio distributions among Root, Leaf and Flower for each RNA type (`mRNA`, `rRNA`, `tRNA`), showing whether modification ratios of all detected WT positive sites differ among tissues within each RNA class.

Re-implementation of `archive/20260902/paper/WT_sites_overview/tissue_ratio` on the current data and the shared figure style (`Figures/shared/plot_style.R` + `plot_tokens.R`).

## Inputs
- `polyARNA/STAT/fisher/mRNA/WT-{Root,leaf,Flower}_sites.positive.tsv` (+ `.Fisher.tsv`)
- `polyARNA/STAT/fisher/rRNA/WT-{Root,leaf,Flower}_sites.positive.tsv` (+ `.Fisher.tsv`)
- `tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.tsv` (+ `.Fisher.tsv`)

Note the sample-name separators: polyARNA uses hyphens (`WT-leaf`), tRNA uses underscores (`WT_Leaf`).

For tRNA, only records with chromosome names matching `^tRNA-` are used (nuclear tRNA; excludes sn/snoRNA mapped to the same reference).

Per site, the treatment depth/gap come from the Fisher table (`Treat_avg_depth`, `Treat_avg_gap`), falling back to the `depth`/`ratio * depth` of the positive table when the join misses.

## Plot
- X axis: RNA type (labels only — no X axis line or tick marks).
- Offset axes: the Y axis line is visually separated from the panel (`guide_prism_offset` + `init_theme_offset`).
- Y axis spans 0–100%; the panel ends exactly at the top tick (data maximum is 100%).
- Within each RNA type, Root, Leaf and Flower are shown as separate bars.
- Bar height: mean modification ratio of all WT positive sites in that tissue and RNA type.
- Error bar: SEM of site-level modification ratios.
- Dots: individual WT positive sites; each dot is one site's modification ratio.
- Significance brackets are drawn above the panel in the top margin (`clip = "off"`), since the axis ends at 100%.
- Colors: shared pastel tissue palette (`PAL_TISSUE`); legend sits on the right with plural display forms (Roots / Leaves / Flowers).
- The dot layer does not generate a legend; the retained legend is generated from the bar layer.
- Canvas: 8.4 x 6.4 in.

## Statistical Test
For each RNA type, pairwise tissue comparisons are performed:
- Root vs Leaf
- Root vs Flower
- Leaf vs Flower

The test input is the site-level modification ratio percentage (`ratio_pct = ratio * 100`). Each pair uses a Welch two-sample t-test (`stats::t.test(x, y, var.equal = FALSE)`), appropriate because the tissues can have different numbers of detected sites and unequal variance.

No multiple-testing correction is applied. Significance labels use the raw Welch t-test p value:
- `***`: p <= 0.001
- `**`: p <= 0.01
- `*`: p <= 0.05
- `n.s.`: p > 0.05

## Outputs
- `tissue_ratio_distribution.pdf`
- `tissue_ratio_distribution.png`
- `tissue_ratio_distribution.svg`
- `tissue_ratio_site_ratios.tsv` — site-level ratios (analysis input)
- `tissue_ratio_summary.tsv` — per RNA type x tissue site counts and pooled ratios
- `tissue_ratio_ttest_tests.tsv` — pairwise Welch t-test results

## Run
```bash
./data.sh
```
