# mRNA_motif — Sequence logos for WT mRNA Fisher positive sites

Sequence logos for wild-type mRNA positive sites across ±2, ±5, and ±10 nt windows (bits scale).

## Input criteria

- Input sites: `polyARNA/STAT/fisher/mRNA/WT-{Root,leaf,Flower}_sites.positive.tsv` (Root: 235 / Leaf: 199 / Flower: 163)
- Sequence extraction: Because Fisher positive tables lack sequence context columns, sequences for all three windows are sliced using 1-based transcript spliced coordinates from Araport11 representative transcript cDNAs (FASTA; sense orientation, no reverse complement). Filtering criteria: uppercase, U→T conversion, complete windows only, no ambiguous bases (N).
- All sites retained without filtering; center base is 100% T across all three window sizes (QC summary in `mrna_motif_center_base_stats.tsv` within each window directory).

## Run

```bash
bash run.sh
```

## Directory layout

```
mRNA_motif/
├── run.sh                        # Master entrypoint
├── extract_window_sequences.py   # Window sequence extraction (FASTA slicing)
├── plot_mrna_motif_seqlogo.R     # Four-panel ggseqlogo (Roots/Leaves/Flowers/All), bits scale
└── +-{2,5,10}/                   # Per-window: sequence table + center base QC + seqlogo in three formats
```

## Key observations

- Position 0 (center base): 100% T (2.0 bits).
- Position -1 exhibits modest T enrichment (Leaf ~0.3, Flower ~0.28, Root ~0.15 bits).
- Position +2 exhibits weak G enrichment (~0.15–0.2 bits, consistent across all three tissues).
- Across ±5 and ±10 nt windows, upstream and downstream signals decay to background, with sequence preferences concentrated within ±2 nt.
