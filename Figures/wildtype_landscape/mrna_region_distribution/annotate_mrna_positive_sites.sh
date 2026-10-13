#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
readonly OUTPUT_DIR="${SCRIPT_DIR}/site_analyse"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly GPPY="/home/HL/software/gppy/gppy/gtf.py"
readonly GTF="${PROJECT_ROOT}/Araport11_fixed.gtf"
readonly REF_GENOME="/home/haoli/reference/Ath/TAIR/TAIR10_chr_all.fasta"
readonly ANNOTATION_SCRIPT="/home/haoli/software/Stript/BASAL_PARISE/bed_annoat.R"

log() {
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*" >&2
}

error() {
  log "ERROR: $*"
  exit 1
}

require_command() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || error "Required command not found: $cmd"
}

require_file() {
  local path="$1"
  [[ -f "$path" ]] || error "Required file not found: $path"
}

main() {
  require_file "$PYTHON_BIN"
  require_command Rscript
  require_command awk
  require_command paste

  require_file "$GPPY"
  require_file "$GTF"
  require_file "$REF_GENOME"
  require_file "$ANNOTATION_SCRIPT"

  mkdir -p "$OUTPUT_DIR"

  shopt -s nullglob
  local input_files=("${PROJECT_ROOT}/polyARNA/STAT/fisher/mRNA/WT"*_sites.positive.tsv)
  shopt -u nullglob

  [[ ${#input_files[@]} -gt 0 ]] || error "No WT mRNA positive.tsv files found under ${PROJECT_ROOT}/polyARNA/STAT/fisher/mRNA"

  for input_file in "${input_files[@]}"; do
    local base_name
    local tsv_file
    local bed_file
    local annotated_bed_file
    local annotated_tsv_file

    base_name="$(basename "$input_file" .tsv)"
    tsv_file="${OUTPUT_DIR}/${base_name}.tsv"
    bed_file="${OUTPUT_DIR}/${base_name}.bed"
    annotated_bed_file="${OUTPUT_DIR}/${base_name}.annoated.bed"
    annotated_tsv_file="${OUTPUT_DIR}/${base_name}.annoated.tsv"

    log "Processing ${base_name}"

    awk 'NR > 1' "$input_file" > "$tsv_file"

    "$PYTHON_BIN" "$GPPY" t2g -g "$GTF" -i "$tsv_file" | \
      awk 'BEGIN{FS="\t"; OFS="\t"} $(NF-2) != "NA" {print $(NF-2), $NF - 1, $NF, ".", ".", $(NF-1)}' \
      > "$bed_file"

    Rscript "$ANNOTATION_SCRIPT" \
      -b "$bed_file" \
      -g "$GTF" \
      -o "$annotated_bed_file"

    printf 'chr\tpos\tstrand\tannoation\n' > "$annotated_tsv_file"
    paste -d $'\t' \
      <(awk -F'\t' '{print $1 "\t" $2 "\t" $3}' "$tsv_file") \
      <(awk -F'\t' '{print $4}' "$annotated_bed_file") \
      >> "$annotated_tsv_file"

    # 中间文件(tsv/bed/annoated.bed)可由本脚本完全再生, 只保留最终注释表
    rm -f "$tsv_file" "$bed_file" "$annotated_bed_file"

    log "Wrote ${annotated_tsv_file}"
  done

  log "All mRNA positive site annotations completed. Outputs are in ${OUTPUT_DIR}"
}

main "$@"
