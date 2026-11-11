#!/usr/bin/env bash
# extract_windows.sh — Extract ±25/±50/±100 nt window sequences around trua1a
# writer-dependent Ψ sites
# Scope:
#   1) case: the unique physical sites of site_info/trua1a_dependent_sites.tsv
#      (transcript/pos/strand in columns 3/4/5);
#   2) all dependent sites per tissue: the per-tissue dependent.tsv of the upstream
#      dependent_betabinomial calling (chr/pos/strand in columns 1/2/3).
# Window logic matches secondary_structure/data_prep.sh exactly: pos is a 1-based
# transcript coordinate, start = pos-w-1 clamped to 0, end = pos+w; bedtools slop -b 0
# clips to the transcript bounds and getfasta -s returns the sense strand.
# Validation: bed/fa row counts must equal the site count; the center base of each ±25
# window (the (w+1)-th base when unclamped) must be T.
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"
readonly DEPENDENT_DIR="${ROOT_DIR}/polyARNA/STAT/dependent_betabinomial/mRNA"
readonly REF_FA="/home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly OUT_DIR="${SCRIPT_DIR}/windows"
readonly -a TISSUES=("Root" "leaf" "Flower")
readonly -a HALF_WINDOWS=(25 50 100)

require_nonempty_file() {
  local file="$1"
  [[ -s "${file}" ]] || {
    echo "[extract_windows] Missing or empty input file: ${file}" >&2
    exit 1
  }
}

# extract_windows <sites_tsv> <col_chr> <col_pos> <col_strand> <prefix>
extract_windows() {
  local sites_tsv="$1" col_chr="$2" col_pos="$3" col_strand="$4" prefix="$5"
  local site_count=$(( $(wc -l < "${sites_tsv}") - 1 ))

  for w in "${HALF_WINDOWS[@]}"; do
    local bed_file="${OUT_DIR}/${prefix}_${w}.bed"
    local fa_file="${OUT_DIR}/${prefix}_${w}.fa"

    awk -F '\t' -v OFS='\t' -v w="${w}" \
        -v col_chr="${col_chr}" -v col_pos="${col_pos}" -v col_strand="${col_strand}" '
      NR > 1 {
        start = $col_pos - w - 1
        if (start < 0) start = 0
        print $col_chr, start, $col_pos + w, $col_chr "_" $col_pos "_" $col_strand, 0, $col_strand
      }
    ' "${sites_tsv}" |
      bedtools slop -i - -g "${REF_FA}.fai" -b 0 > "${bed_file}"

    bedtools getfasta -fi "${REF_FA}" -bed "${bed_file}" -s -bedOut |
      awk -F '\t' '{ print ">" $4; print toupper($7) }' > "${fa_file}"

    local bed_count fa_count
    bed_count="$(wc -l < "${bed_file}")"
    fa_count=$(( $(wc -l < "${fa_file}") / 2 ))
    if [[ "${bed_count}" -ne "${site_count}" || "${fa_count}" -ne "${site_count}" ]]; then
      echo "[extract_windows] Row count mismatch for ${prefix} at w=${w}: bed=${bed_count} fa=${fa_count} sites=${site_count}" >&2
      exit 1
    fi
    echo "[extract_windows] ${prefix} ±${w} nt: ${bed_count} windows -> ${fa_file}"
  done

  # Center-base check for ±25 windows (clipped windows locate the center by pos-start0)
  "${PYTHON_BIN}" - "${OUT_DIR}/${prefix}_25.bed" "${OUT_DIR}/${prefix}_25.fa" <<'PYEOF'
import sys
bed = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
seqs = {}
name = None
for line in open(sys.argv[2]):
    line = line.rstrip("\n")
    if line.startswith(">"):
        name = line[1:]
    else:
        seqs[name] = line
bad = 0
for chrom, start, end, name, _score, strand in bed:
    tx, pos, _ = name.rsplit("_", 2)
    center = int(pos) - int(start) - 1  # 0-based center index within the window
    if seqs[name][center] != "T":
        bad += 1
        print(f"[extract_windows] center base not T: {name} -> {seqs[name][center]}",
              file=sys.stderr)
if bad:
    sys.exit(f"[extract_windows] {bad} windows failed center-base check in {sys.argv[2]}")
PYEOF
}

require_nonempty_file "${REF_FA}"
require_nonempty_file "${REF_FA}.fai"
mkdir -p "${OUT_DIR}"

# 1) case sites
case_tsv="${SCRIPT_DIR}/site_info/trua1a_dependent_sites.tsv"
require_nonempty_file "${case_tsv}"
extract_windows "${case_tsv}" 3 4 5 "case_windows"

# 2) all dependent sites per tissue
for tissue in "${TISSUES[@]}"; do
  dep_tsv="${DEPENDENT_DIR}/${tissue}/trua1a/trua1a-${tissue}.dependent.tsv"
  require_nonempty_file "${dep_tsv}"
  extract_windows "${dep_tsv}" 1 2 3 "${tissue}_dependent_windows"
done

echo "[extract_windows] Done."
