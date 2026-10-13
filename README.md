# PUS Writer Dependent sites Analysis

This repository contains the analysis code used to process BID-seq transcriptome data for the PUS writer-dependent pseudouridylation study. The main mapping and candidate-site calling workflow is implemented in [`BID-transcript.script/workflows/polyA_mapping.smk`](BID-transcript.script/workflows/polyA_mapping.smk).

The workflow is designed for paired Mock/Treat BID libraries. It preprocesses FASTQ files, maps reads to transcript references, corrects local gap signals, filters candidate T sites, and applies replicate-aware Fisher testing to identify positive sites.

## Repository layout

- [`BID-transcript.script/`](BID-transcript.script): Snakemake workflows and helper programs for BID transcript mapping and site calling; see its own [`README.md`](BID-transcript.script/README.md) for the full layout.
- [`BID-transcript.script/workflows/`](BID-transcript.script/workflows): Snakemake workflows — `polyA_mapping.smk` (main poly(A) RNA BID-seq workflow), `totalRNA_mapping.smk`, and `polyA_mapping-adjustDeletion.smk`.
- [`BID-transcript.script/bin/`](BID-transcript.script/bin): helper programs invoked by the workflows (the config's `script_dir`), including `Pseudo_callsites.R` for Fisher testing and positive-site calling.
- [`BID-transcript.script/downstream/`](BID-transcript.script/downstream): standalone downstream tools, e.g. `call_dependent.py` for PUS-dependent site calling and `fit_calibration.py` for ratio calibration.
- [`BID-transcript.script/config.example.yaml`](BID-transcript.script/config.example.yaml): example configuration file for running the workflow.
- `paper/`: manuscript-related notes.

## Upstream code

The helper programs [`adjustGap`](BID-transcript.script/bin/adjustGap), [`realignGap`](BID-transcript.script/bin/realignGap), [`samFilter`](BID-transcript.script/bin/samFilter), and [`cpup`](BID-transcript.script/bin/cpup) are from the [y9c/pseudoU-BIDseq](https://github.com/y9c/pseudoU-BIDseq) project and are used here as part of the BID-seq mapping and gap-processing workflow.

## Pipeline overview

For each replicate-level FASTQ pair, `polyA_mapping.smk` runs:

1. adapter and length filtering with `trim_galore`;
2. duplicate removal with `fastp --dedup`;
3. fixed UMI/adapter-region clipping with `fastp`;
4. primary mapping to the first configured reference with `bowtie2`;
5. optional secondary mapping of unmapped reads with `STAR` when two reference types are configured;
6. forward-strand BAM filtering, gap realignment, sorting, `calmd`, and mismatch filtering;
7. `samtools mpileup` and conversion to per-site count tables;
8. combined gap adjustment across all replicates for each reference type;
9. T-site filtering with local sequence context extraction;
10. Mock-vs-Treat Fisher testing with `Pseudo_callsites.R`.

The reference types are read from `pipeline.ref` in YAML order. The first reference is used for the Bowtie2 stage. If a second reference is provided, the workflow treats it as a STAR stage and requires `paths.index_star`.

## Requirements

The commands called by the workflow must be available in the runtime environment:

- Snakemake
- `trim_galore`
- `fastp`
- `bowtie2`
- `STAR`
- `samtools`
- `bedtools`
- `csvtk`
- Python
- R

The helper binaries/scripts in `BID-transcript.script/bin` are also required: `samFilter`, `realignGap`, `cpup`, `adjustGap`, `combine_mpileup_tsv.py`, and `Pseudo_callsites.R`.

## Configure a run

Start from the example file:

```bash
cp BID-transcript.script/config.example.yaml config.yaml
```

Edit `config.yaml` for your dataset:

- `paths.script_dir`: absolute path to `BID-transcript.script/bin`.
- `paths.pseudo_callsites_r`: absolute path to `BID-transcript.script/bin/Pseudo_callsites.R`.
- `paths.index_bowtie2`: Bowtie2 index basename for the first reference type.
- `paths.index_star`: STAR genome directory; required only when two reference types are configured.
- `env.python`: Python executable used for `combine_mpileup_tsv.py`.
- `env.r`: R executable used to run `Pseudo_callsites.R`.
- `pipeline.ref`: one or two named reference FASTA paths. Reference FASTA files must have `.fai` indexes for the T-site context step.
- `pipeline.t_site_windows`: number of bases on each side of the candidate T site to extract as context.
- `samples`: biological sample definitions. Each sample must contain matched `mock` and `treat` replicates. Replicate names must match between groups, and every `file_index` must be unique.

Relative FASTQ paths in `samples.*.*.*.read1/read2` are resolved relative to the Snakemake working directory passed by `--directory`.

## Run the workflow

The analysis run command follows the reference script at `/home/HL/project/shuangli/PUS/PUScaller/Mapping_pipeline/BID-mod/polyARNA/run.sh`. A portable version of that launcher is:

```bash
#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

snakemake \
  -s "/path/to/BID-transcript.script/workflows/polyA_mapping.smk" \
  --configfile "$SCRIPT_DIR/config.yaml" \
  --directory "$SCRIPT_DIR" \
  --cores 200 \
  --rerun-incomplete \
  "$@"
```

In practice, create a run directory that contains `config.yaml` and `raw_data/`, then run:

```bash
snakemake \
  -s /path/to/BID-transcript.script/workflows/polyA_mapping.smk \
  --configfile config.yaml \
  --directory . \
  --cores 200 \
  --rerun-incomplete
```

Additional Snakemake options can be appended as needed, for example `-n` for a dry run or `--printshellcmds` for command logging.

## Main outputs

The workflow writes outputs under the run directory:

- `fix.fastq/`: trimmed and deduplicated FASTQ files.
- `umi_extract.fastq/`: reads after fixed UMI/adapter-region clipping.
- `BAM/`: mapped, realigned, filtered BAM files and mapping reports.
- `STAT/mpileup/<ref_type>/`: replicate-level mpileup TSV files and `COMBINED-<ref_type>.tsv`.
- `STAT/adjustGap/<ref_type>/`: gap-adjusted tables and filtered T-site inputs.
- `STAT/fisher/<ref_type>/`: `<sample>_sites.Fisher.tsv` and `<sample>_sites.positive.tsv`.
- `logs/`: per-replicate read-count and mapping summary tables.

The positive-site tables in `STAT/fisher/<ref_type>/` are the main candidate pseudouridylation site calls for downstream analysis.
