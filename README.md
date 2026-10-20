# PUS writer-dependent pseudouridylation — analysis code and data

<p align="center">
  <a href="README.md">English</a> · <a href="README.zh-CN.md">中文</a>
</p>

This repository accompanies the manuscript:

> **[Title TBD]**
> [Authors TBD]. *Journal TBD*. Article link: TBD.

It contains the analysis code and processed data used in the paper:

- the Snakemake workflows that map BID-seq libraries and call pseudouridylation
  (Ψ) sites,
- the downstream scripts that identify PUS writer-dependent sites and
  temperature-responsive sites,
- the exact configuration of every analysis dataset,
- the per-figure plotting code together with its input tables.

## Repository layout

| Path | Content |
|---|---|
| [`BID-transcript.script/`](BID-transcript.script/) | Mapping and site-calling workflows (`workflows/`), helper programs (`bin/`), and downstream callers (`downstream/`). See its own [README](BID-transcript.script/README.md) and [CONFIG_GUIDE](BID-transcript.script/CONFIG_GUIDE.md) for full documentation. |
| [`configs/`](configs/) | The five analysis configurations used in this study (one YAML per dataset; see **Datasets and configurations**). |
| [`Figures/`](Figures/) | Figure-by-figure plotting code. Each leaf directory holds a `data_prep.sh`, the committed input tables (`*.tsv`/`*.bed`), a `plot.R`, and the rendered PDF/PNG/SVG, so every plot can be re-rendered from the committed tables without re-running the pipelines. Shared styling lives in `Figures/shared/`. |

Raw FASTQ files and large intermediates (BAMs, mpileup tables) are not part of
the repository; see **Data availability**.

## Analysis pipeline

BID-seq (bisulfite-induced deletion sequencing) libraries are sequenced as
Mock/Treat pairs: in Treat libraries, bisulfite treatment converts Ψ into a
Ψ–bisulfite adduct that induces deletion (gap) signatures at Ψ sites during
reverse transcription, while unmodified uridines give near-zero background;
Mock libraries are the untreated controls
([Li et al., *Nat Plants* 2025](https://doi.org/10.1038/s41477-024-01894-7)).

Each dataset is processed end-to-end by one of two Snakemake workflows —
[`polyA_mapping.smk`](BID-transcript.script/workflows/polyA_mapping.smk) for
poly(A)+ RNA libraries and
[`totalRNA_mapping.smk`](BID-transcript.script/workflows/totalRNA_mapping.smk)
for tRNA-enriched / total-RNA libraries (both workflows analyze only R2) —
through the following stages:

1. **Preprocessing** — adapter, quality and length trimming with `trim_galore`;
   PCR duplicate removal (`fastp --dedup` for poly(A) libraries, `seqkit rmdup`
   for total-RNA libraries); fixed UMI/adapter-region clipping with `fastp`,
   parameterized per read end (`umi_extract` in the config).
2. **Mapping** — `bowtie2` against the first configured reference (an rRNA
   reference for poly(A) libraries, a unique tRNA/sn/snoRNA reference for
   tRNA-enriched libraries). When two references are configured, unmapped reads
   are subsequently aligned with `STAR` to the Araport11 representative cDNA
   models.
3. **Signal extraction** — the mapped BAMs are processed into per-site count
   tables: `samFilter` strand-aware filtering → `realignGap` local gap
   realignment → `samtools calmd` mismatch annotation and filtering →
   `samtools mpileup` per-site depth and gap counts.
4. **Gap correction** — `adjustGap` jointly re-estimates gap signals across all
   replicates of a dataset, parallelized over reference chunks
   (`adjust_gap.chunk_bp` / `max_sequences`).
5. **Site calling** — candidate T (Ψ) sites are filtered by depth
   (`filter_t_sites`), then `bin/Pseudo_callsites.R` runs replicate-aware
   Mock-vs-Treat Fisher tests with per-replicate pass criteria (the `filter`
   block in the config), producing `<sample>_sites.Fisher.tsv` and
   `<sample>_sites.positive.tsv` under `STAT/fisher/<ref_type>/`.
6. **Downstream calls** — PUS writer-dependent sites are called from the
   positive-site tables with
   [`downstream/call_dependent_methylKit.R`](BID-transcript.script/downstream/call_dependent_methylKit.R);
   temperature-responsive sites with
   [`downstream/call_temperature_response_methylKit.R`](BID-transcript.script/downstream/call_temperature_response_methylKit.R).
   In parallel, the workflows emit raw per-gene expression matrices under
   `count/<ref_type>/` (Mock-group and all-samples).

The helper programs `adjustGap`, `realignGap`, `samFilter` and `cpup` originate
from the [y9c/pseudoU-BIDseq](https://github.com/y9c/pseudoU-BIDseq) project.

## Site-calling methods

The three statistical callers invoked at stages 5 and 6 above share one
representation of the data: for every candidate T (Ψ) site and every
replicate, the gap-corrected pileup provides the sequencing depth (reads
spanning the site) and the gap count (reads carrying a deletion at the
site); the ratio of gap count to depth estimates that replicate's Ψ
fraction. This is formally
identical to a bisulfite-seq measurement — coverage plus methylated-read
counts — which is why the two differential callers below reuse the methylKit
differential-methylation framework
([Akalin et al., *Genome Biol* 2012](https://doi.org/10.1186/gb-2012-13-9-r87)).

### Identification of Ψ sites — [`bin/Pseudo_callsites.R`](BID-transcript.script/bin/Pseudo_callsites.R)

**Merging and pre-filtering.** For every Mock/Treat sample pair, the
per-replicate site tables are merged on site coordinates (chr, position,
strand); sites absent from a replicate enter the merge with zero counts, and
sites whose group-mean depth falls below 10 in either group are discarded.

**Group-level test.** Bisulfite-induced deletion enrichment is tested with a
one-sided Fisher's exact test on the replicate-pooled 2×2 table of gap versus
non-gap reads (Treat against Mock, alternative "greater"), and p-values are
Benjamini–Hochberg-corrected across sites (retained at adjusted p < 0.01).

**Per-replicate criteria.** In addition, each replicate must individually
support the call:

- a Treat gap count ≥ 5;
- a Treat gap ratio ≥ 5%;
- a Treat gap ratio at least twice that of the paired Mock replicate.
  Replicates are compared pairwise, Mock replicate *i* against Treat
  replicate *i*; when the two groups have unequal replicate numbers, each
  Treat replicate is instead compared against the Mock group mean.

**Positive call.** A site is called positive when it passes the group-level
test and at least `min_pass_reps` replicates pass the per-replicate criteria
— 3 of 3 for the three-replicate datasets, 2 of 2 for the two-replicate
chloroplast dataset. In the positive-site table, the reported ratio and depth
are averaged over the passing replicates only.

### Identification of PUS writer-dependent sites — [`downstream/call_dependent_methylKit.R`](BID-transcript.script/downstream/call_dependent_methylKit.R)

**The idea.** If a Ψ site is installed by a particular PUS enzyme, mutating
that enzyme should make the modification disappear — the Treat-library gap
ratio at the site should fall back towards the Mock background. The script
therefore asks one question per site — *is the modification fraction
significantly lower in the mutant than in the WT?* — and answers it
separately for each of the two independent mutant alleles of the writer gene.

**Why a DNA-methylation tool can test this.** Per site and per replicate,
BID-seq produces just two numbers: how many reads cover the site (depth) and
how many of them carry a deletion (gap). That is exactly what bisulfite-seq
produces for a cytosine (coverage and methylated reads), so the standard
differential-methylation package methylKit can be reused as-is by reading
the gap ratio as the "methylation frequency".

**What the script does, step by step, per mutant allele:**

1. **Candidates.** Only the sites called positive in the WT dataset are
   examined — in the mutants, we only ask whether WT-modified sites are
   still modified.
2. **Depth gate.** Every replicate (WT, mutant Mock and mutant Treat) must
   cover the site with depth ≥ 10; one failing replicate drops the site.
   This stops "not sequenced" from masquerading as "not modified" (the
   adjusted tables record uncovered sites as depth = 0) and keeps only
   complete cases. Mock libraries witness coverage only; the signal compared
   is the Treat gap ratio.
3. **Test.** WT and mutant Treat replicates go through methylKit
   (`methRead` → `unite` → `calculateDiffMeth`): a coverage-weighted
   logistic regression of the gap ratio on genotype (WT = 1, mutant = 0),
   with overdispersion estimated as in McCullagh and Nelder — the genotype
   effect is referred to an F distribution when φ > 1 (chi-square otherwise)
   — and p-values adjusted across sites with Benjamini–Hochberg (BH). The
   BH q-value (two-sided FDR) is used for the call, with direction taken
   from meth.diff: q < 0.05 and meth.diff > 0 (higher modification in WT)
   are both required; the opposite direction never passes.
4. **Effect size.** Δ = WT_Ratio − Mut_Ratio, each term the arithmetic mean
   of the replicate Treat gap ratios. It is computed from the input tables,
   independently of methylKit, so the reported difference always matches
   the raw ratios.
5. **Call.** Dependent in one allele: BH q < 0.05, meth.diff > 0, **and**
   Δ ≥ 0.05.
   The final writer-dependent set keeps only sites dependent in **both**
   alleles, so a site cannot pass on the peculiarities of a single line;
   the mutant ratio reported there averages the two lines.

Outputs: per-allele `<allele>.dependent-sites.tsv`, the combined
`<gene>.dependent.tsv`, and a summary JSON.

### Identification of temperature-responsive sites — [`downstream/call_temperature_response_methylKit.R`](BID-transcript.script/downstream/call_temperature_response_methylKit.R)

Temperature-responsive sites are called with the same methylKit workflow,
comparing the same WT material grown at the two temperatures (27 °C vs
21 °C).

**Candidates and depth gate.** Candidates are the union of the positive
sites called at either temperature (a site positive at only one temperature
is precisely a candidate for a temperature-dependent change in modification),
merged upstream of the call. The same per-replicate depth gate (≥ 10 in
every Mock and Treat replicate of both temperatures) is applied — this also
excludes sites whose transcript is not expressed at one temperature, which
would otherwise masquerade as a modification difference.

**Test and call.** Because a temperature response can go in either
direction, the two-sided test is used without a direction filter: a site is
temperature-responsive when the BH q-value (two-sided FDR) < 0.05 and
|Δ| ≥ 0.05 (group ratios are arithmetic means of the replicate Treat gap
ratios), and the direction (increase/decrease at 27 °C) is recorded per site.

Outputs: `<group2>_vs_<group1>.all-sites.tsv`, `.response.tsv`, and a summary
JSON with up/down counts.

## Datasets and configurations

[`configs/`](configs/) holds the exact YAML configuration of each of the five
analysis datasets in this study (in the working tree, each file is symlinked as
`config.yaml` into its pipeline run directory):

| Config | Workflow | Dataset |
|---|---|---|
| [`polyARNA.yaml`](configs/polyARNA.yaml) | polyA_mapping | poly(A)+ RNA from WT and 14 mutant lines of seven TRUA/PUS genes, each in root, leaf and flower |
| [`tRNA.yaml`](configs/tRNA.yaml) | totalRNA_mapping | tRNA-enriched libraries, same genotypes and tissues |
| [`Temperature.yaml`](configs/Temperature.yaml) | polyA_mapping | poly(A)+ RNA of WT and two trua1b lines grown at 21 °C and 27 °C |
| [`Temperature_tRNA.yaml`](configs/Temperature_tRNA.yaml) | totalRNA_mapping | tRNA-enriched libraries, same temperature design |
| [`Chloroplast.yaml`](configs/Chloroplast.yaml) | totalRNA_mapping | chloroplast RNA |

### Configuration parameters

All five configs share one schema:

- **`paths`** — absolute locations of the helper programs (`script_dir`,
  `pseudo_callsites_r`) and of the aligner indexes (`index_bowtie2` for the
  first reference, `index_star` required only with two references).
- **`env`** — the Python and R executables used by workflow scripts.
- **`pipeline.ref`** — one or two named reference FASTAs, in mapping order: the
  first is the bowtie2 stage, the (optional) second is the STAR stage. The same
  names key the output directories (`STAT/<ref_type>/`, `count/<ref_type>/`).
- **`pipeline.t_site_windows`** — bases of sequence context extracted on each
  side of a candidate T site.
- **`pipeline.expression`** — optional raw expression matrices per reference
  type; `count_mode` is `all_alignments` (alignment-level, for multi-mapping
  RNAs such as rRNA/tRNA) or `primary` (read-level, for mRNA).
- **`pipeline.adjust_gap`** — chunking of the joint gap-correction step:
  `chunk_bp` per chunk, batching small transcripts up to `max_sequences`.
- **`pipeline.filter`** — site-calling thresholds:
  - `min_depth`: per-replicate depth for the candidate pool;
  - `min_avg_depth`: group-level mean depth;
  - `max_padj`: BH-adjusted Fisher p-value cutoff;
  - `min_gap`: Treat gap count per passing replicate;
  - `min_treat_ratio`: Treat gap ratio per passing replicate;
  - `ratio_fold`: Treat/Mock ratio fold-change;
  - `min_pass_reps`: how many replicates must pass — 3 for all
    three-replicate datasets, 2 for the two-replicate Chloroplast dataset.
- **`read_layout` / `rmdup` / `trim_galore_length`** — paired vs single-end
  input, whether to deduplicate, and the minimum read length kept after
  trimming.
- **`pipeline.umi_extract`** — fixed-length clipping of read termini with
  `fastp` (`trim_front`/`trim_tail` = bases removed from the 5′/3′ end). The
  clipping lengths follow from the library structure: every read carries a
  fixed 8 nt UMI plus the `TATAGGG` sequence, and the first few low-quality
  5′ bases are removed as well.
  - **Poly(A) libraries**: 6 such bases are removed, so 21 nt are clipped
    from the 5′ end (8 + 7 + 6) and 6 nt from the 3′ end, and the paired-end
    libraries are clipped complementarily as `R1: 6/21, R2: 21/6`.
  - **tRNA-enriched / total-RNA libraries**: reads are short, so only the
    first 2 low-quality bases are removed; the 5′ end of R2 is
    `8 bp UMI + TATAGGG + 2 bp`, hence `R2: 17/0` (mirrored as `R1: 0/17`).
  - For both library types only R2 enters mapping and downstream analysis.
- **`samples`** — per-sample replicate definitions:
  `samples.<sample>.{mock,treat}.<rep>.{file_index,read1,read2}`, with FASTQ
  paths relative to the run directory's `raw_data/`. An optional per-replicate
  **`stages`** key restricts a replicate to a branch: `stages: [count]`
  replicates only contribute to mapping/expression matrices and never enter the
  Fisher site-calling branch, so adding them does not invalidate the joint
  gap-correction results (the default, `[count, fisher]`, enters both).

### Running a dataset

Each run directory holds its `config.yaml` (symlinked into `configs/`) and a
`raw_data/` directory of FASTQ files, and is executed with:

```bash
snakemake \
  -s /path/to/BID-transcript.script/workflows/polyA_mapping.smk \
  --configfile config.yaml \
  --directory . \
  --cores 200 \
  --rerun-incomplete
```

(use `totalRNA_mapping.smk` for the total-RNA datasets). Main outputs are
`STAT/mpileup/`, `STAT/adjustGap/`, `STAT/fisher/` (positive-site tables),
`count/` (expression matrices) and `logs/` (per-replicate library statistics).

## Requirements

### Mapping and site-calling workflows

- **Orchestration**: [Snakemake](https://snakemake.github.io).
- **Command-line tools**: `trim_galore` (which internally calls cutadapt and
  FastQC), `fastp`, `seqkit`, `bowtie2`, `STAR`, `samtools`, `bedtools`,
  `csvtk`.
- **Python 3**: [`polars`](https://pola.rs/),
  [`pysam`](https://pysam.readthedocs.io/), PyYAML.
- The helper programs `adjustGap`, `realignGap`, `samFilter` and `cpup` ship
  with the repository under `BID-transcript.script/bin/` and need no separate
  installation.

### Downstream differential calls

- The Bioconductor package
  [`methylKit`](https://bioconductor.org/packages/methylKit/) (PUS
  writer-dependent and temperature-responsive site calling).

## Data availability

The raw sequencing data generated in this study have been deposited under
accession **[TBD]** (link TBD). Processed per-site tables underlying every
figure are committed alongside the plotting code under [`Figures/`](Figures/).

## Citation

If you use this code or data, please cite:

> [Citation TBD]
