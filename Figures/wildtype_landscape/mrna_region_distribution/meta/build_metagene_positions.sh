#!/usr/bin/env bash

# Build the site coordinate table for the mRNA metagene density figure.
#
# WT positive sites from the three tissues (polyARNA/STAT/fisher/mRNA/
# WT-*_sites.positive.tsv; coordinates are 1-based spliced positions on the
# representative transcripts) are projected onto a 5'UTR-CDS-3'UTR metagene:
# each transcript's three regions are rescaled to fixed fractional segments
# (15% / 70% / 15% by default); a site's x = region start + relative position
# within the segment x region fraction. Region boundaries come directly from
# utr5_len/cds_len/tx_len of gppy txinfo_basic:
#   5'UTR = [1, utr5_len], CDS = (utr5_len, utr5_len+cds_len], 3'UTR = rest.
# Sites on non-coding transcripts (cds_len = 0), absent from the GTF, or out
# of range are dropped and counted; they do not enter the density.
#
# Outputs: meta/mrna_metagene_sites.tsv
#          (Tissue/transcript/pos/strand/ratio/depth/Region/Region_len/x)
#          meta/mrna_metagene_sites_stats.tsv (kept/dropped counts per tissue)

set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# meta sits four levels below the repo root:
# meta -> mrna_region_distribution -> wildtype_landscape -> Figures -> repo
readonly PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
readonly OUTPUT_DIR="${SCRIPT_DIR}"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly GPPY="/home/HL/software/gppy/gppy/gtf.py"
readonly GTF="${PROJECT_ROOT}/Araport11_fixed.gtf"
readonly FISHER_DIR="${PROJECT_ROOT}/polyARNA/STAT/fisher/mRNA"
readonly OUT_TSV="${OUTPUT_DIR}/mrna_metagene_sites.tsv"
readonly OUT_STATS="${OUTPUT_DIR}/mrna_metagene_sites_stats.tsv"

# Fixed fractional lengths of the three metagene segments:
# 5'UTR 0-15%, CDS 15-85%, 3'UTR 85-100%
readonly U5_FRAC=0.15
readonly CDS_FRAC=0.70

# Tissue spelling in file prefixes -> display name in the output tables
readonly TISSUES=("Root" "leaf" "Flower")
readonly TISSUE_LABELS=("Root" "Leaf" "Flower")

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

error() {
  log "ERROR: $*"
  exit 1
}

main() {
  [[ -f "$PYTHON_BIN" ]] || error "Required file not found: $PYTHON_BIN"
  [[ -f "$GPPY" ]] || error "Required file not found: $GPPY"
  [[ -f "$GTF" ]] || error "Required file not found: $GTF"

  # Deliberately not local: the EXIT trap fires only after main returns,
  # when local variables are already out of scope
  tmp_dir="$(mktemp -d "${OUTPUT_DIR}/.metagene.XXXXXX")"
  trap 'rm -rf "$tmp_dir"' EXIT

  # Per-transcript region boundaries: tx_name tx_len cds_len utr5_len
  # (non-coding transcripts included; used for the drop classification)
  "$PYTHON_BIN" "$GPPY" txinfo_basic -g "$GTF" 2>/dev/null \
    | awk -F'\t' 'BEGIN{OFS="\t"} FNR == 1 {next} {print $1, $6, $7, $8}' \
    > "$tmp_dir/tx_regions.tsv"
  log "Transcript region map: $(wc -l < "$tmp_dir/tx_regions.tsv") transcripts"

  printf 'Tissue\ttranscript\tpos\tstrand\tratio\tdepth\tRegion\tRegion_len\tx\n' > "$OUT_TSV"
  printf 'Tissue\tInput\tMapped\tDrop_tx_absent\tDrop_noncoding\tDrop_out_of_range\n' > "$OUT_STATS"

  for i in "${!TISSUES[@]}"; do
    tissue="${TISSUES[$i]}"
    label="${TISSUE_LABELS[$i]}"
    input_file="${FISHER_DIR}/WT-${tissue}_sites.positive.tsv"
    [[ -f "$input_file" ]] || error "Required file not found: $input_file"

    log "Processing WT-${tissue} (${label})"

    # Site classification and metagene coordinates; the count row is written
    # to a separate file and appended to the stats table afterwards
    awk -F'\t' -v OFS='\t' \
      -v tissue="$label" -v u5f="$U5_FRAC" -v cdsf="$CDS_FRAC" \
      -v r5="5'UTR" -v r3="3'UTR" \
      -v statsfile="$tmp_dir/stats.${label}.tsv" '
      NR == FNR { txlen[$1] = $2; cds[$1] = $3; u5[$1] = $4; next }
      {
        n_input++
        t = $1; p = $2 + 0
        if (!(t in txlen)) { n_notx++; next }
        if (cds[t] == 0) { n_noncoding++; next }
        c = u5[t] + cds[t]
        if (p <= u5[t]) {
          r = r5; off = p - 1; len = u5[t]; start = 0; frac = u5f
        } else if (p <= c) {
          r = "CDS"; off = p - u5[t] - 1; len = cds[t]; start = u5f; frac = cdsf
        } else if (p <= txlen[t]) {
          r = r3; off = p - c - 1; len = txlen[t] - c; start = u5f + cdsf; frac = 1 - u5f - cdsf
        } else {
          n_oor++; next
        }
        if (len <= 0) { n_oor++; next }
        printf "%s\t%s\t%d\t%s\t%s\t%s\t%s\t%d\t%.6f\n", \
          tissue, t, p, $3, $4, $5, r, len, start + off / len * frac
        n_mapped++
      }
      END {
        print tissue, n_input, n_mapped, n_notx, n_noncoding, n_oor > statsfile
      }' "$tmp_dir/tx_regions.tsv" "$input_file" >> "$OUT_TSV"

    cat "$tmp_dir/stats.${label}.tsv" >> "$OUT_STATS"
    log "WT-${tissue}: $(awk -F'\t' '{print $3 " mapped / " $2 " input"}' "$tmp_dir/stats.${label}.tsv")"
  done

  log "Wrote ${OUT_TSV}"
  log "Wrote ${OUT_STATS}"
}

main "$@"
