#!/usr/bin/env bash
set -euo pipefail

IFS=$'\n\t'

# Sequence logos (y axis in bits) of Fisher-positive WT mRNA Ψ sites for
# three windows: ±2 / ±5 / ±10 nt. The Fisher-positive tables lack a context
# column, so all three windows are sliced by transcript coordinates from
# Araport11 representative-transcript cDNA (FASTA). Each window's sequence
# table, center-base QC, and logo figure go into its own +-{2,5,10}/
# subfolder.
#
# Usage: bash run.sh
# Outputs: +-{2,5,10}/mrna_motif_sequences.tsv
#          +-{2,5,10}/mrna_motif_center_base_stats.tsv
#          +-{2,5,10}/mrna_motif_seqlogo.{pdf,png,svg}

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# mRNA_motif -> wildtype_landscape -> Figures -> repo root
readonly PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
readonly FISHER_MRNA_DIR="${PROJECT_ROOT}/polyARNA/STAT/fisher/mRNA"
readonly REF_FASTA="/home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly RSCRIPT_BIN="/home/HL/miniforge3/envs/analysis/bin/Rscript"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

main() {
  [[ -f "$PYTHON_BIN" ]] || { echo "Missing python: $PYTHON_BIN" >&2; exit 1; }
  [[ -f "$RSCRIPT_BIN" ]] || { echo "Missing Rscript: $RSCRIPT_BIN" >&2; exit 1; }
  [[ -f "$REF_FASTA" ]] || { echo "Missing FASTA: $REF_FASTA" >&2; exit 1; }
  [[ -d "$FISHER_MRNA_DIR" ]] || { echo "Missing dir: $FISHER_MRNA_DIR" >&2; exit 1; }

  log "Extracting window sequences (±2 / ±5 / ±10)"
  "$PYTHON_BIN" "$SCRIPT_DIR/extract_window_sequences.py" \
    "$FISHER_MRNA_DIR" "$REF_FASTA" "$SCRIPT_DIR"

  for window in 2 5 10; do
    local width=$((2 * window + 1))
    local window_dir="${SCRIPT_DIR}/+-${window}"
    log "Drawing seqlogo for ±${window} (${width}-mer)"
    "$RSCRIPT_BIN" "$SCRIPT_DIR/plot_mrna_motif_seqlogo.R" \
      "${window_dir}/mrna_motif_sequences.tsv" "$window_dir" "$width"
  done

  log "Done"
  log "Figures: $SCRIPT_DIR/+-{2,5,10}/mrna_motif_seqlogo.{pdf,png,svg}"
}

main "$@"
