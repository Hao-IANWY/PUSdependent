#!/usr/bin/env bash
# data_prep.sh — extract TRUA-dependent tRNA 5-nt contexts by genotype × tissue
# 职责：生成 plot.R 的唯一数据接口 trna_motif_sequences.tsv 和计数表。
# 只保留 chr 以 tRNA- 开头的行（参考序列混有 sn/snoRNA，必须剔除）。
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
readonly DEPENDENT_DIR="${ROOT_DIR}/tRNA/STAT/dependent_betabinomial/ncRNA"
readonly GENOTYPES=(trua1a trua1b trua2 trua3 trua4 trua5 trua6)
readonly TISSUES=(Root leaf Flower)
readonly SEQUENCE_TSV="${SCRIPT_DIR}/trna_motif_sequences.tsv"
readonly COUNT_TSV="${SCRIPT_DIR}/trna_motif_sequence_counts.tsv"
readonly TMP_SEQUENCE_TSV="${SCRIPT_DIR}/.trna_motif_sequences.tsv.tmp"

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

append_one_combo() {
  local genotype="$1"
  local tissue="$2"
  local input_file="${DEPENDENT_DIR}/${tissue}/${genotype}/${genotype}-${tissue}.dependent.tsv"

  require_nonempty_file "${input_file}"

  awk -F '\t' -v OFS='\t' \
      -v genotype="${genotype}" -v tissue="${tissue}" '
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
    $(col["chr"]) !~ /^tRNA-/ { next }
    {
      context = toupper($(col["context"]))
      gsub(/[^ACGTU]/, "", context)
      if (length(context) != 5) next

      sequence = context
      gsub(/U/, "T", sequence)
      printf("%s\t%s\t%s\t%s\t%s\t%s\n", genotype, tissue,
             $(col["chr"]), $(col["pos"]), $(col["strand"]), sequence)
    }
  ' "${input_file}" >> "${TMP_SEQUENCE_TSV}"
}

main() {
  printf 'genotype\ttissue\tchr\tpos\tstrand\tsequence\n' > "${TMP_SEQUENCE_TSV}"

  local genotype tissue
  for genotype in "${GENOTYPES[@]}"; do
    for tissue in "${TISSUES[@]}"; do
      append_one_combo "${genotype}" "${tissue}"
    done
  done

  require_nonempty_file "${TMP_SEQUENCE_TSV}"
  mv "${TMP_SEQUENCE_TSV}" "${SEQUENCE_TSV}"

  awk -F '\t' -v OFS='\t' '
    NR == 1 {
      print "genotype", "tissue", "Sequence_Count"
      next
    }
    { count[$1 SUBSEP $2]++ }
    END {
      n_g = split("trua1a trua1b trua2 trua3 trua4 trua5 trua6", gs, " ")
      n_t = split("Root leaf Flower", ts, " ")
      for (gi = 1; gi <= n_g; gi++) {
        for (ti = 1; ti <= n_t; ti++) {
          key = gs[gi] SUBSEP ts[ti]
          print gs[gi], ts[ti], ((key in count) ? count[key] : 0)
        }
      }
    }
  ' "${SEQUENCE_TSV}" > "${COUNT_TSV}"

  local row_count
  row_count="$(( $(wc -l < "${SEQUENCE_TSV}") - 1 ))"
  echo "[data_prep] Generated: ${SEQUENCE_TSV} (${row_count} sequences)"
  echo "[data_prep] Generated: ${COUNT_TSV}"
}

main "$@"
