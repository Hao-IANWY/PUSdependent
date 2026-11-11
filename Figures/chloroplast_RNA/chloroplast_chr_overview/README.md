# chloroplast_chr_overview — Circular chloroplast genome (ChrC) Ψ site overview

## Purpose

Generate a circular genome plot of the Arabidopsis chloroplast chromosome based on `Araport11_fixed.gtf` `ChrC` annotations: the outer ring shows chloroplast rRNA / tRNA / mRNA gene positions, and the inner ring shows WT leaf Ψ modification fractions as radial bars.

## Inputs

- `Araport11_fixed.gtf`
- `polyARNA/STAT/fisher/rRNA/WT-leaf_sites.positive.tsv`
- `tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.positive.tsv`
- `polyARNA/STAT/fisher/mRNA/WT-leaf_sites.positive.tsv`

## Coordinate Mapping

- mRNA sites: `ATCG*.1` transcript IDs are matched directly to `ChrC` transcript annotations in the GTF and converted to `ChrC` genomic coordinates via cumulative exon offsets (reversed on the minus strand).
- rRNA (`Ch*`) and tRNA (`tRNA-Pt-*`) site tables use transcript-internal coordinates. `plot.R` contains an internal mapping table from `Ch*` / `tRNA-Pt-*` to `ATCG*` gene IDs, converted using GTF strand and exon structures.
- The output table `chloroplast_chr_overview_sites.tsv` retains fields including `source_chr`, `gene_id`, `pos`, `chr_pos`, and `ratio_pct` to verify coordinate mapping.

## Data Summary

- ChrC features: 88 mRNA, 8 rRNA, 37 tRNA
- WT leaf sites mapped to ChrC: 4 rRNA, 63 tRNA, 18 mRNA
- 1 unmapped tRNA site: `tRNA-Pt-Ser-TGA-1:74` (falls outside the exon range of this gene), written to `chloroplast_chr_overview_unmapped_sites.tsv` with a warning

## Outputs

- `chloroplast_chr_overview.pdf/png/svg`: Circular genome plot
- `chloroplast_chr_overview_annotated.pdf/png/svg`: Circular genome plot with mRNA gene ID annotations
- `chloroplast_chr_overview_features.tsv`: ChrC features used for plotting
- `chloroplast_chr_overview_sites.tsv`: WT leaf sites mapped to ChrC
- `chloroplast_chr_overview_unmapped_sites.tsv`: Unmapped sites
- `chloroplast_chr_overview_feature_summary.tsv`: Summary of feature counts
- `chloroplast_chr_overview_site_summary.tsv`: Summary of site counts and mean modification fractions

Legend: Three distinct colors represent tRNA / mRNA / rRNA on both outer ring genes and inner ring sites; radial bar length is proportional to modification fraction, with outer black dots marking site positions.

## Usage

```bash
bash data.sh
```
