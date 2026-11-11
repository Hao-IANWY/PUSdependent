#!/usr/bin/env bash

# Compute the total spliced length of each mRNA-region class in the reference
# (representative gene models); these lengths are the denominator of the
# per-1000-nt modification-density figure.
#
# The mRNA reference is the Araport11 representative gene model set (one
# transcript per gene), so isoforms cannot be double counted. Region lengths
# come from gppy txinfo (spliced lengths derived from exon structure,
# utr3 = tx_len - cds_len - utr5_len) and sum(tx_len) is cross-checked
# against the summed length column of the reference .fai.
#
# Outputs: mrna_region_distribution_region_lengths.tsv
#          (Region/Length_nt; one row each for 5'UTR/CDS/3'UTR)

set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly GPPY="/home/HL/software/gppy/gppy/gtf.py"
readonly GTF="${PROJECT_ROOT}/Araport11_fixed.gtf"
readonly MRNA_REF="/home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model"
readonly OUTPUT_TSV="${SCRIPT_DIR}/mrna_region_distribution_region_lengths.tsv"

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
  [[ -f "${MRNA_REF}.fai" ]] || error "Required file not found: ${MRNA_REF}.fai"

  # Deliberately not local: the EXIT trap fires only after main returns, when
  # local variables are already out of scope
  tmp_dir="$(mktemp -d "${SCRIPT_DIR}/.region_lengths.XXXXXX")"
  trap 'rm -rf "$tmp_dir"' EXIT

  # Reference transcript id set + total spliced length (.fai column 2)
  cut -f1 "${MRNA_REF}.fai" > "$tmp_dir/ref_ids.txt"
  local n_ref fai_total
  n_ref="$(wc -l < "$tmp_dir/ref_ids.txt")"
  fai_total="$(awk '{s += $2} END {print s}' "${MRNA_REF}.fai")"
  log "Reference transcripts: ${n_ref}, total spliced length (.fai): ${fai_total}"

  # Per-transcript spliced region lengths:
  # tx_name gene_id chrom strand nexon tx_len cds_len utr5_len utr3_len
  # (txinfo_basic is used because txinfo raises KeyError on organellar genes
  # while parsing attributes)
  "$PYTHON_BIN" "$GPPY" txinfo_basic -g "$GTF" > "$tmp_dir/tx_info.tsv"

  # Restrict to the reference set and sum per region; non-coding transcripts
  # with utr3_len == 0 naturally contribute no CDS length
  awk -F'\t' '
    NR==FNR { ref[$1]=1; next }
    FNR==1 { next }
    ($1 in ref) { tx+=$6; cds+=$7; u5+=$8; u3+=$9; n++ }
    END {
      # QC row (matched_tx/sum_tx_len) goes to stderr; region-length rows
      # go to stdout
      printf "%d\t%d\t%d\t%d\n", n, tx, cds, u5+u3 > "/dev/stderr"
      printf "5\047UTR\t%d\nCDS\t%d\n3\047UTR\t%d\n", u5, cds, u3
    }' "$tmp_dir/ref_ids.txt" "$tmp_dir/tx_info.tsv" \
    > "$tmp_dir/region_lengths.tsv" 2> "$tmp_dir/sums.txt"

  local matched sum_tx
  matched="$(cut -f1 "$tmp_dir/sums.txt")"
  sum_tx="$(cut -f2 "$tmp_dir/sums.txt")"
  log "Matched transcripts in GTF: ${matched} / ${n_ref}"
  log "sum(tx_len) = ${sum_tx}  vs  .fai total = ${fai_total}"
  [[ "$sum_tx" == "$fai_total" ]] || log "WARNING: tx_len sum differs from .fai total"
  [[ "$matched" == "$n_ref" ]] || log "WARNING: $((n_ref - matched)) reference transcript(s) missing from GTF"

  # List reference ids absent from the GTF (e.g. organellar transcripts) for
  # inspection
  awk -F'\t' 'NR==FNR { if (FNR > 1) seen[$1]=1; next } !($1 in seen) { print }' \
    "$tmp_dir/tx_info.tsv" "$tmp_dir/ref_ids.txt" > "$tmp_dir/unmatched_ids.txt"
  if [[ -s "$tmp_dir/unmatched_ids.txt" ]]; then
    log "Unmatched reference id(s): $(paste -sd, "$tmp_dir/unmatched_ids.txt")"
  fi

  printf 'Region\tLength_nt\n' > "$OUTPUT_TSV"
  cat "$tmp_dir/region_lengths.tsv" >> "$OUTPUT_TSV"

  log "Wrote ${OUTPUT_TSV}"
}

main "$@"
