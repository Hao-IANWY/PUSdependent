#!/usr/bin/env bash

# Build the per-tissue expressed pool of assayable uridines for WT mRNA and
# annotate it by mRNA region.
#
# The expressed pool = every uridine actually covered by sequencing in that
# tissue's samples and therefore testable for modification; it is the
# denominator of per-region modification rates. Background: the tables
# polyARNA/STAT/fisher/mRNA/WT-*_sites.Fisher.tsv do not contain every tested
# site -- Pseudo_callsites.R already filters pooled Fisher p.adj < 0.01 before
# writing them. The pool is therefore rebuilt from the per-replicate
# candidate files using the same entry criteria as that script:
#   1. STAT/adjustGap/mRNA/*.filter_T.tsv: sites with reference base T and
#      depth >= 10 in a given replicate;
#   2. union across all 6 replicates (Mock x3 + Treat x3), with depth set to
#      0 for replicates missing a site;
#   3. summed Mock depth >= 30 and summed Treat depth >= 30 (i.e. group mean
#      >= min_avg_depth = 10).
# Pool sites then pass through exactly the same annotation path as the
# positive sites (gppy t2g to genomic coordinates + bed_annoat.R region
# annotation), keeping numerator (positive sites) and denominator (expressed
# pool) on identical region calling criteria.
#
# Outputs: site_analyse/WT-<Tissue>.expressed_pool.annoated.tsv (chr/pos/strand/annoation)

set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
readonly OUTPUT_DIR="${SCRIPT_DIR}/site_analyse"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly GPPY="/home/HL/software/gppy/gppy/gtf.py"
readonly GTF="${PROJECT_ROOT}/Araport11_fixed.gtf"
readonly ANNOTATION_SCRIPT="/home/haoli/software/Stript/BASAL_PARISE/bed_annoat.R"
readonly FILTER_T_DIR="${PROJECT_ROOT}/polyARNA/STAT/adjustGap/mRNA"
readonly MIN_SUM_DEPTH=30
readonly T2G_CHUNKS=16

readonly TISSUES=("Root" "leaf" "Flower")

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

error() {
  log "ERROR: $*"
  exit 1
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || error "Required file not found: $path"
}

main() {
  require_file "$PYTHON_BIN"
  require_file "$GPPY"
  require_file "$GTF"
  require_file "$ANNOTATION_SCRIPT"

  mkdir -p "$OUTPUT_DIR"

  for tissue in "${TISSUES[@]}"; do
    local prefix="WT-${tissue}"
    local rep_files=(
      "${FILTER_T_DIR}/Mock-${prefix}.rep1.mRNA.filter_T.tsv"
      "${FILTER_T_DIR}/Mock-${prefix}.rep2.mRNA.filter_T.tsv"
      "${FILTER_T_DIR}/Mock-${prefix}.rep3.mRNA.filter_T.tsv"
      "${FILTER_T_DIR}/Treat-${prefix}.rep1.mRNA.filter_T.tsv"
      "${FILTER_T_DIR}/Treat-${prefix}.rep2.mRNA.filter_T.tsv"
      "${FILTER_T_DIR}/Treat-${prefix}.rep3.mRNA.filter_T.tsv"
    )
    for f in "${rep_files[@]}"; do
      require_file "$f"
    done

    local pool_tsv="${OUTPUT_DIR}/${prefix}.expressed_pool.tsv"
    local pool_bed="${OUTPUT_DIR}/${prefix}.expressed_pool.bed"
    local annotated_bed="${OUTPUT_DIR}/${prefix}.expressed_pool.annoated.bed"
    local annotated_tsv="${OUTPUT_DIR}/${prefix}.expressed_pool.annoated.tsv"

    log "Building expressed pool for ${prefix}"

    # Union across replicates with group-mean depth >= 10 (Mock/Treat, 3
    # replicates each; a missing replicate counts 0, so the summed threshold
    # = 3 x 10)
    awk -F'\t' -v min_sum="$MIN_SUM_DEPTH" 'BEGIN{OFS="\t"}
      FNR == 1 {next}
      {
        key = $1 "|" $2
        if (!(key in str)) { str[key] = $3 }
        if (ARGIND <= 3) { mock[key] += $6 } else { treat[key] += $6 }
      }
      END {
        for (key in str) {
          if (mock[key] >= min_sum && treat[key] >= min_sum) {
            split(key, a, "|")
            print a[1], a[2], str[key]
          }
        }
      }' "${rep_files[@]}" > "$pool_tsv"

    log "${prefix} pool size: $(wc -l < "$pool_tsv")"

    # Append the row number as an id in column 4 (t2g reads only the first
    # two columns, chr/pos): rows whose t2g mapping fails (transcript absent
    # from the GTF) are dropped consistently on both the pool and bed sides
    # by id, keeping rows aligned with the annotation output
    local numbered_tsv="${OUTPUT_DIR}/${prefix}.expressed_pool.numbered.tsv"
    awk 'BEGIN{OFS="\t"}{print $0, NR}' "$pool_tsv" > "$numbered_tsv"

    # t2g is slow on millions of rows; run chunks in parallel and concatenate
    # chunk outputs in chunk order to preserve the row order of pool.tsv
    local chunk_dir
    chunk_dir="$(mktemp -d "${OUTPUT_DIR}/.expressed_pool_chunks.XXXXXX")"
    split -n "l/${T2G_CHUNKS}" -d "$numbered_tsv" "${chunk_dir}/chunk."

    local pids=()
    for chunk in "${chunk_dir}"/chunk.*; do
      "$PYTHON_BIN" "$GPPY" t2g -g "$GTF" -i "$chunk" 2>/dev/null | \
        awk -v chunk="$chunk" 'BEGIN{FS="\t"; OFS="\t"}
          $(NF-2) != "NA" {
            print $(NF-2), $NF - 1, $NF, ".", ".", $(NF-1) > (chunk ".bed");
            print $4 > (chunk ".ids");
            next
          }
          { print $1 > (chunk ".dropped") }' &
      pids+=($!)
    done
    for pid in "${pids[@]}"; do
      wait "$pid" || error "t2g chunk failed for ${prefix}"
    done

    local mapped_tsv="${OUTPUT_DIR}/${prefix}.expressed_pool.mapped.tsv"
    cat "${chunk_dir}"/chunk.*.ids | sort -n > "${chunk_dir}/keep.ids"
    awk 'NR==FNR { keep[$1]=1; next } keep[$4]' "${chunk_dir}/keep.ids" "$numbered_tsv" > "$mapped_tsv"
    cat "${chunk_dir}"/chunk.*.bed > "$pool_bed"
    if compgen -G "${chunk_dir}/chunk.*.dropped" > /dev/null; then
      sort -u "${chunk_dir}"/chunk.*.dropped > "${OUTPUT_DIR}/${prefix}.expressed_pool.unmapped_transcripts.txt"
      log "${prefix}: dropped $(wc -l < "${OUTPUT_DIR}/${prefix}.expressed_pool.unmapped_transcripts.txt") transcript(s) absent from GTF"
    fi
    rm -rf "$chunk_dir" "$numbered_tsv"

    local n_mapped n_bed
    n_mapped="$(wc -l < "$mapped_tsv")"
    n_bed="$(wc -l < "$pool_bed")"
    [[ "$n_mapped" == "$n_bed" ]] || error "${prefix}: row mismatch between mapped pool (${n_mapped}) and pool.bed (${n_bed})"

    Rscript "$ANNOTATION_SCRIPT" \
      -b "$pool_bed" \
      -g "$GTF" \
      -o "$annotated_bed"

    printf 'chr\tpos\tstrand\tannoation\n' > "$annotated_tsv"
    paste -d $'\t' \
      <(awk -F'\t' '{print $1 "\t" $2 "\t" $3}' "$mapped_tsv") \
      <(awk -F'\t' '{print $4}' "$annotated_bed") \
      >> "$annotated_tsv"

    # Intermediates (pool/mapped/bed/annoated.bed) are fully regenerable by
    # this script; keep only the final annotation table
    rm -f "$pool_tsv" "$mapped_tsv" "$pool_bed" "$annotated_bed"

    log "Wrote ${annotated_tsv}"
  done

  log "All tested pools built and annotated. Outputs are in ${OUTPUT_DIR}"
}

main "$@"
