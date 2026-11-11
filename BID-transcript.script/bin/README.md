# bin/

Utility scripts for per-site signal processing and candidate Ψ-site calling
from BID-seq mpileup data.

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
