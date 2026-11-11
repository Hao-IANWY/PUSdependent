#!/usr/bin/env bash
# data_prep.sh — Fold the four window/site groups with ViennaRNA and build plotting TSVs
# Inputs:
#   windows/case_windows_{25,50,100}.{bed,fa}: unique physical TRUA1A-dependent sites
#   windows/{Root,Leaf,Flower}_dependent_windows_{25,50,100}.{bed,fa}: per-tissue sites
# Outputs:
#   <group>_structure_metrics.tsv and <group>_structure_unpaired_profile.tsv
#   combined trua1a_groups_structure_metrics.tsv / trua1a_groups_structure_unpaired_profile.tsv
# Model: ViennaRNA 2.7.0 (RNA.fold_compound), GU wobble allowed, no pseudoknots, 21 °C.
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"
readonly DEPENDENT_DIR="${ROOT_DIR}/polyARNA/STAT/dependent_betabinomial/mRNA"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly FOLD_SCRIPT="${SCRIPT_DIR}/../secondary_structure/fold_structure.py"
readonly SITE_DIR="${SCRIPT_DIR}/site_sets"
readonly WINDOWS_DIR="${SCRIPT_DIR}/windows"
readonly TEMPS="21"

require_nonempty_file() {
  local file="$1"
  [[ -s "${file}" ]] || {
    echo "[data_prep] Missing or empty input file: ${file}" >&2
    exit 1
  }
}

require_nonempty_file "${PYTHON_BIN}"
require_nonempty_file "${FOLD_SCRIPT}"
mkdir -p "${SITE_DIR}"

# Build the minimal site tables expected by the validated folding script. The context
# column is retained so fold_structure.py can verify coordinates, strand, and FASTA identity.
printf 'genotype\ttissue\ttranscript\tpos\tstrand\tcontext\tset\n' > "${SITE_DIR}/unique_sites.tsv"
awk -F '\t' -v OFS='\t' 'NR > 1 {print "trua1a", "Union", $3, $4, $5, $6, "unique"}' \
  "${SCRIPT_DIR}/site_info/trua1a_dependent_sites.tsv" >> "${SITE_DIR}/unique_sites.tsv"

for tissue in Root leaf Flower; do
  display_tissue="$(printf '%s' "${tissue:0:1}" | tr '[:lower:]' '[:upper:]')${tissue:1}"
  input="${DEPENDENT_DIR}/${tissue}/trua1a/trua1a-${tissue}.dependent.tsv"
  require_nonempty_file "${input}"
  out="${SITE_DIR}/${display_tissue}_sites.tsv"
  printf 'genotype\ttissue\ttranscript\tpos\tstrand\tcontext\tset\n' > "${out}"
  awk -F '\t' -v OFS='\t' -v tissue="${display_tissue}" '
    NR > 1 {print "trua1a", tissue, $1, $2, $3, $4, tissue}
  ' "${input}" >> "${out}"
done

# Fold each group. The sibling fold_structure.py performs a context check at every site
# and extracts base-neighborhood unpaired probabilities from the partition function (BPP).
fold_group() {
  local group="$1" sites="$2" prefix="$3"
  require_nonempty_file "${SITE_DIR}/${sites}"
  "${PYTHON_BIN}" "${FOLD_SCRIPT}" "${SCRIPT_DIR}" \
    --sites "site_sets/${sites}" \
    --bed-prefix "windows/${prefix}" \
    --out-prefix "${group}_structure" \
    --temps "${TEMPS}"
}

fold_group unique unique_sites.tsv case_windows
fold_group Root Root_sites.tsv Root_dependent_windows
fold_group Leaf Leaf_sites.tsv leaf_dependent_windows
fold_group Flower Flower_sites.tsv Flower_dependent_windows

# Combine the four groups without changing any fold-derived values.
combine_tsv() {
  local suffix="$1" out="$2" first=1
  for group in unique Root Leaf Flower; do
    local file="${SCRIPT_DIR}/${group}_structure_${suffix}"
    require_nonempty_file "${file}"
    if (( first )); then
      cat "${file}" > "${out}"
      first=0
    else
      tail -n +2 "${file}" >> "${out}"
    fi
  done
}
combine_tsv metrics.tsv "${SCRIPT_DIR}/trua1a_groups_structure_metrics.tsv"
combine_tsv unpaired_profile.tsv "${SCRIPT_DIR}/trua1a_groups_structure_unpaired_profile.tsv"

profile_rows=$(( $(wc -l < "${SCRIPT_DIR}/trua1a_groups_structure_unpaired_profile.tsv") - 1 ))
metric_rows=$(( $(wc -l < "${SCRIPT_DIR}/trua1a_groups_structure_metrics.tsv") - 1 ))
echo "[data_prep] Wrote: ${SCRIPT_DIR}/trua1a_groups_structure_metrics.tsv (${metric_rows} rows)"
echo "[data_prep] Wrote: ${SCRIPT_DIR}/trua1a_groups_structure_unpaired_profile.tsv (${profile_rows} rows)"
