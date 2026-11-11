# PUS writer-dependent pseudouridylation: analysis code and data

<p align="center">
  <a href="README.md">English</a> · <a href="README.zh-CN.md">Chinese</a>
</p>

This repository contains the analysis code and processed data associated with the following manuscript:

> **[Title TBD]**
> [Authors TBD]. *Journal TBD*. Article link: TBD.

This repository includes:

- **Snakemake workflows**: Automated pipelines for mapping BID-seq libraries and identifying pseudouridine (Ψ) sites.
- **Downstream analysis scripts**: Statistical methods for identifying PUS writer-dependent and temperature-responsive Ψ sites.
- **Analysis configurations**: Comprehensive YAML configuration files for all datasets used in this study.
- **Visualization code**: Scripts and input tables for generating figures, enabling reproducible rendering of all plots.

## Repository layout

| Path | Content |
|---|---|
| [`BID-transcript.script/`](BID-transcript.script/) | Core analysis pipelines (`workflows/`), auxiliary programs (`bin/`), and downstream statistical scripts (`downstream/`). Refer to the [README](BID-transcript.script/README.md) and [CONFIG_GUIDE](BID-transcript.script/CONFIG_GUIDE.md) for detailed documentation. |
| [`configs/`](configs/) | Dataset-specific configuration files (one YAML per dataset). See **Datasets and configurations** for details. |
| [`Figures/`](Figures/) | Figure-specific plotting scripts and source data. Each subdirectory contains a `data_prep.sh`, the required input tables (`*.tsv`/`*.bed`), a `plot.R` script, and rendered outputs (PDF/PNG/SVG). |

*Note: Raw FASTQ files and large intermediate files (e.g., BAMs, mpileup tables) are excluded; see **Data availability**.*

## Analysis pipeline

BID-seq (bisulfite-induced deletion sequencing) libraries are analyzed in Mock/Treat pairs. In Treat libraries, bisulfite treatment converts Ψ into a Ψ–bisulfite adduct, which induces site-specific deletions during reverse transcription, whereas untreated Mock libraries provide minimal background ([Li et al., *Nat Plants* 2025](https://doi.org/10.1038/s41477-024-01894-7)).

Each dataset is processed end-to-end using one of two primary Snakemake workflows:
- [`polyA_mapping.smk`](BID-transcript.script/workflows/polyA_mapping.smk) for poly(A)+ RNA libraries.
- [`totalRNA_mapping.smk`](BID-transcript.script/workflows/totalRNA_mapping.smk) for tRNA-enriched / total-RNA libraries.

Both workflows analyze only R2 reads through the following stages:

1. **Preprocessing**: Adapter trimming and quality filtering with `trim_galore`, followed by PCR duplicate removal (`fastp --dedup` for poly(A), `seqkit rmdup` for total-RNA). Fixed UMI/adapter-region clipping is performed via `fastp`, with parameters defined per read end in the configuration (`umi_extract`).
2. **Mapping**: Reads are aligned using `bowtie2` against the primary reference (rRNA for poly(A) libraries; tRNA/sn/snoRNA for tRNA-enriched libraries). If a secondary reference is configured, unmapped reads are further aligned using `STAR` against the Araport11 representative cDNA models.
3. **Signal extraction**: BAM files are processed into site-specific count tables: `samFilter` (strand-aware filtering) → `realignGap` (local gap realignment) → `samtools calmd` (mismatch annotation and filtering) → `samtools mpileup` (per-site depth and gap quantification).
4. **Gap correction**: `adjustGap` re-estimates gap signals jointly across all replicates, parallelized over reference chunks (`adjust_gap.chunk_bp` / `max_sequences`).
5. **Site calling**: Candidate Ψ sites are filtered by depth (`filter_t_sites`). Subsequently, `bin/Pseudo_callsites.R` performs replicate-aware Mock-vs-Treat Fisher tests based on criteria defined in the config, outputting results to `STAT/fisher/<ref_type>/` (e.g., `<sample>_sites.Fisher.tsv` and `<sample>_sites.positive.tsv`).
6. **Downstream analysis**: PUS writer-dependent sites and temperature-responsive sites are identified using [`downstream/call_dependent_betabinomial.R`](BID-transcript.script/downstream/call_dependent_betabinomial.R) and [`downstream/call_temperature_response_betabinomial.R`](BID-transcript.script/downstream/call_temperature_response_betabinomial.R), respectively. Simultaneously, the workflows generate raw per-gene expression matrices in `count/<ref_type>/`.

The helper programs `adjustGap`, `realignGap`, `samFilter`, and `cpup` are provided by the [y9c/pseudoU-BIDseq](https://github.com/y9c/pseudoU-BIDseq) project.

## Statistical methodology

The statistical callers share a unified data representation: for each candidate Ψ site, the gap-corrected pileup provides total depth (reads spanning the site) and gap counts (reads with a deletion). The ratio of gap counts to total depth serves as the estimated Ψ fraction, modeled via a beta-binomial generalized linear mixed model (GLMM) to account for replicate-level overdispersion.

### Identification of Ψ sites ([`bin/Pseudo_callsites.R`](BID-transcript.script/bin/Pseudo_callsites.R))

**Merging and pre-filtering**: Sample pairs (Mock/Treat) are merged by coordinate (chr, position, strand). Sites with group-mean depth < 10 are excluded.

**Group-level testing**: Bisulfite-induced deletion enrichment is tested using a one-sided Fisher's exact test (Treat vs. Mock; alternative: "greater"). P-values are Benjamini–Hochberg (BH) corrected across sites; sites with adjusted p < 0.01 are retained.

**Per-replicate criteria**: A site is called positive if it passes the group-level test and satisfies the following criteria in at least `min_pass_reps` replicates:
- Treat gap count ≥ 5;
- Treat gap ratio ≥ 5%;
- Treat gap ratio ≥ 2× that of the paired Mock replicate.
*(For datasets with unequal replicate numbers, Treat replicates are compared against the Mock group mean.)*

The reported ratio and depth are calculated as the mean across qualifying replicates.

### PUS writer-dependent sites ([`downstream/call_dependent_betabinomial.R`](BID-transcript.script/downstream/call_dependent_betabinomial.R))

**Objective**: To determine if the modification fraction significantly decreases in the mutant compared to the Wild Type (WT).

**Workflow**:
1. **Candidates**: Only sites identified as positive in the WT are analyzed.
2. **Depth gate**: Requires depth ≥ 10 in all WT and mutant (Mock/Treat) replicates.
3. **Statistical model**: A beta-binomial GLMM is fitted using `glmmTMB`. First, a site-level random intercept estimates global dispersion (θ), which is then fixed. Each site is tested for genotype effect (`cbind(gap, depth − gap) ~ genotype`) via a likelihood-ratio χ²(1) test.
4. **Significance**: One-sided p-values (based on the β sign) are BH-corrected; sites with q < 0.05 and Δ (WT_Ratio − Mut_Ratio) ≥ the effect-size floor (0.05 for mRNA and tRNA, 0.10 for rRNA) are considered dependent.
5. **Final call**: Only sites dependent in **both** independent mutant alleles are retained.

### Temperature-responsive sites ([`downstream/call_temperature_response_betabinomial.R`](BID-transcript.script/downstream/call_temperature_response_betabinomial.R))

Temperature-responsive sites are identified using the same beta-binomial framework, comparing WT samples grown at 21 °C vs. 27 °C. Candidates include the union of positive sites at either temperature. A two-sided likelihood-ratio test is applied: sites with BH q < 0.05 and |Δ| ≥ 0.05 are classified as temperature-responsive.

## Datasets and configurations

[`configs/`](configs/) contains the YAML configurations for the five datasets:

| Config | Workflow | Description |
|---|---|---|
| [`polyARNA.yaml`](configs/polyARNA.yaml) | polyA_mapping | poly(A)+ RNA: WT and 14 mutant lines of 7 PUS genes (root, leaf, flower) |
| [`tRNA.yaml`](configs/tRNA.yaml) | totalRNA_mapping | tRNA-enriched libraries: same genotypes and tissues |
| [`Temperature.yaml`](configs/Temperature.yaml) | polyA_mapping | poly(A)+ RNA: WT and 2 *trua1b* lines grown at 21 °C and 27 °C |
| [`Temperature_tRNA.yaml`](configs/Temperature_tRNA.yaml) | totalRNA_mapping | tRNA-enriched libraries: same temperature design |
| [`Chloroplast.yaml`](configs/Chloroplast.yaml) | totalRNA_mapping | Chloroplast RNA |

### Configuration parameters
Key configuration fields:
- `paths`: Paths to helper programs and alignment indexes.
- `pipeline.ref`: Reference FASTAs, specified in alignment order.
- `pipeline.filter`: Site-calling thresholds (e.g., `min_depth`, `max_padj`, `min_gap`, `min_pass_reps`).
- `pipeline.umi_extract`: Clipping parameters for read termini, derived from library-specific UMI/adapter structures.
- `samples`: Replicate definitions, FASTQ paths, and optional `stages` (e.g., `[count]`) to restrict processing branches.

## Requirements

### Analysis pipeline
- **Orchestration**: [Snakemake](https://snakemake.github.io).
- **Tools**: `trim_galore`, `fastp`, `seqkit`, `bowtie2`, `STAR`, `samtools`, `bedtools`, `csvtk`.
- **Dependencies**: Python 3 (`polars`, `pysam`, `PyYAML`).

### Downstream analysis
- R packages: `glmmTMB`, `optparse`, `data.table`, `jsonlite`.

## Data availability

Raw sequencing data will be deposited under accession **[TBD]**. Processed site-specific tables are available within [`Figures/`](Figures/).

## Citation
If you use this code or data, please cite:
> [Citation TBD]
