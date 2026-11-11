# BID-transcript Analysis Pipeline

Scripts for BID-seq (bisulfite-induced deletion sequencing) analysis of
pseudouridylation in transcriptomic and non-coding RNAs. In bisulfite-treated
libraries (Treat), pseudouridine (Ψ) forms a bisulfite adduct that induces
deletion signatures during reverse transcription; paired untreated libraries
(Mock) provide the background. Read 2 is analyzed throughout.

## bin/

Per-site signal processing and candidate Ψ-site calling:

- `Pseudo_callsites.R`: Calls candidate Ψ sites from deletion signals using a replicate-aware one-sided Fisher test.
- `adjustDeletion.py`: Adjusts raw deletion signals for library-specific biases.
- `write_adjust_gap_chunks.py`: Splits mpileup data into chunks for joint gap re-estimation.
- `concatenate_adjusted_chunks.py`: Reassembles gap-adjusted chunks into per-sample tables.
- `combine_mpileup_tsv.py`: Merges per-sample mpileup TSV files into one table.
- `split_mpileup_chunks.awk`: Splits an mpileup TSV into fixed-size chunks by chromosome.
- `summarize_deletion_signals.py`: Summarizes deletion signals across samples.
- `count_feature_alignments.py`: Counts alignments overlapping genomic features.
- `build_library_log.py`: Builds a per-library processing log.
- `summarize_processing_logs.py`: Aggregates per-library processing logs into one report.

Vendored executables (`adjustGap`, `realignGap`, `samFilter`, `cpup`) are
third-party binaries and are kept unmodified.

## downstream/

Statistical calling and annotation of writer-dependent Ψ sites:

- `call_dependent.py`: Wrapper that drives the differential Ψ-site callers.
- `call_dependent_betabinomial.R`: Beta-binomial GLMM (glmmTMB) test for writer-dependent Ψ sites.
- `call_dependent_methylKit.R`: Differential Ψ-site analysis using methylKit.
- `call_temperature_response_betabinomial.R`: Beta-binomial GLMM test for temperature-responsive Ψ sites.
- `convert_to_transcript_coord.py`: Converts genomic coordinates to transcript coordinates.
- `extract_siteinfo.py`: Extracts per-site annotations for called Ψ sites.

### tRNA Annotation

- `add_sprinzl_pos.py`: Maps tRNA site coordinates to canonical Sprinzl positions (0-76).

## workflows/

- `polyA_mapping.smk`: Snakemake workflow for poly(A)-enriched RNA libraries.
- `totalRNA_mapping.smk`: Snakemake workflow for total RNA libraries.
