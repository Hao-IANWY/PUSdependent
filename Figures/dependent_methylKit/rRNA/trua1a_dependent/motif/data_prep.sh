#!/usr/bin/env bash
# data_prep.sh — extract TRUA1A-dependent rRNA 5-nt contexts by tissue
# 职责：生成 plot.R 的唯一数据接口 trua1a_motif_sequences.tsv 和计数表。
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"
readonly DEPENDENT_DIR="${ROOT_DIR}/polyARNA/STAT/dependent_methylKit/rRNA"
readonly GENOTYPE="trua1a"
readonly SEQUENCE_TSV="${SCRIPT_DIR}/trua1a_motif_sequences.tsv"
readonly COUNT_TSV="${SCRIPT_DIR}/trua1a_motif_sequence_counts.tsv"
readonly TMP_SEQUENCE_TSV="${SCRIPT_DIR}/.trua1a_motif_sequences.tsv.tmp"

cleanup() {
  rm -f "${TMP_SEQUENCE_TSV}"
}
trap cleanup EXIT

require_nonempty_file() {
  local file="$1"
  [[ -s "${file}" ]] || {
    echo "[data_prep] Missing or empty input file: ${file}" >&2
    exit 1
  }
}

append_one_tissue() {
  local tissue="$1"
  local input_file="${DEPENDENT_DIR}/${tissue}/${GENOTYPE}/${GENOTYPE}-${tissue}.dependent.tsv"

  require_nonempty_file "${input_file}"

  awk -F '\t' -v OFS='\t' \
      -v genotype="${GENOTYPE}" -v tissue="${tissue}" '
    BEGIN {
      split("chr,pos,strand,context", required_list, ",")
      for (i in required_list) required[required_list[i]] = 1
    }
    FNR == 1 {
      for (i = 1; i <= NF; i++) col[$i] = i
      for (name in required) {
        if (!(name in col)) {
          printf("Required column %s is missing in %s\n", name, FILENAME) > "/dev/stderr"
          exit 1
        }
      }
      next
    }
    {
      context = toupper($(col["context"]))
      gsub(/[^ACGTU]/, "", context)
      if (length(context) != 5) next

      sequence = context
      sub(/U/, "T", sequence)
      gsub(/U/, "T", sequence)
      printf("%s\t%s\t%s\t%s\t%s\t%s\n", genotype, tissue,
             $(col["chr"]), $(col["pos"]), $(col["strand"]), sequence)
    }
  ' "${input_file}" >> "${TMP_SEQUENCE_TSV}"
}

main() {
  printf 'genotype\ttissue\tchr\tpos\tstrand\tsequence\n' > "${TMP_SEQUENCE_TSV}"

  append_one_tissue "Root"
  append_one_tissue "leaf"
  append_one_tissue "Flower"

  require_nonempty_file "${TMP_SEQUENCE_TSV}"
  mv "${TMP_SEQUENCE_TSV}" "${SEQUENCE_TSV}"

  awk -F '\t' -v OFS='\t' '
    NR == 1 {
      print "Tissue", "Sequence_Count"
      next
    }
    { count[$2]++ }
    END {
      print "Root", count["Root"] + 0
      print "leaf", count["leaf"] + 0
      print "Flower", count["Flower"] + 0
    }
  ' "${SEQUENCE_TSV}" > "${COUNT_TSV}"

  local row_count
  row_count="$(( $(wc -l < "${SEQUENCE_TSV}") - 1 ))"
  echo "[data_prep] Generated: ${SEQUENCE_TSV} (${row_count} sequences)"
  echo "[data_prep] Generated: ${COUNT_TSV}"
}

main "$@"
