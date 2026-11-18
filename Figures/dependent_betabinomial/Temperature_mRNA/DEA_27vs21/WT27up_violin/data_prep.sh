#!/usr/bin/env bash
# ============================================================================
# data_prep.sh — WT 27°C 上调基因在 {WT, trua1b_1, trua1b_2} 的 27°C TPM
#
# 基因集：WT_27C_vs_21C_DESeq2.csv 中 padj < 0.05 且 log2FoldChange >= 1
#         的 Upregulated 基因（默认 891 个；与上级火山图阈值一致）
# 表达量：
#   rate = raw_count / transcript_length_kb
#   TPM  = rate / sum(rate across ALL transcripts in the sample) × 1,000,000
#   raw count     : count_27vs21.tsv 中各样品 27°C rep3/4/5 的 primary
#                   mapped-read 计数（来自 realign.filtered.BAM）
#   transcript bp : 比对所用 cDNA 参考 FASTA 的 .fai 序列长度
#   TPM 分母      : 该样品全部转录本的 rate 总和，而非仅上调基因
#
# 输出：
#   WT27up_expr_27C.tsv          长表，含 raw count、长度、TPM 分母和 TPM
#   WT27up_TPM_normalization.tsv 九个样品的 TPM 归一化分母
#
# 用法：
#   bash data_prep.sh
#   LOGFC_FILTER=1 PADJ_FILTER=0.05 bash data_prep.sh   # 覆盖基因集阈值
#   REF_FAI=/path/to/reference.fai bash data_prep.sh     # 覆盖长度索引
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."                       # 上级 DEA_27vs21 目录
COUNT_TSV="${BASE_DIR}/count_27vs21.tsv"
REF_FAI="${REF_FAI:-/home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model.fai}"
LOGFC_FILTER="${LOGFC_FILTER:-1}"
PADJ_FILTER="${PADJ_FILTER:-0.05}"
OUT_TSV="${SCRIPT_DIR}/WT27up_expr_27C.tsv"
NORM_TSV="${SCRIPT_DIR}/WT27up_TPM_normalization.tsv"

[ -f "$COUNT_TSV" ] || { echo "[data_prep] ERROR: count matrix not found: $COUNT_TSV" >&2; exit 1; }
[ -f "$REF_FAI" ] || { echo "[data_prep] ERROR: reference FAI not found: $REF_FAI" >&2; exit 1; }

# ---- WT 上调基因列表（padj < PADJ_FILTER 且 log2FC >= LOGFC_FILTER）----
GENE_LIST="$(mktemp)"
trap 'rm -f "$GENE_LIST"' EXIT
awk -F',' -v lfc="$LOGFC_FILTER" -v padj="$PADJ_FILTER" '
  NR > 1 && $7 != "NA" && $7 != "" && $7+0 < padj && $3+0 >= lfc { print $1 }
' "${BASE_DIR}/WT/WT_27C_vs_21C_DESeq2.csv" | sort -u > "$GENE_LIST"

N_GENES=$(wc -l < "$GENE_LIST")
[ "$N_GENES" -gt 0 ] || { echo "[data_prep] ERROR: no up genes selected" >&2; exit 1; }

# ---- 一次读取长度索引与 count 矩阵，计算 9 个样品的 TPM ----
awk -F'\t' -v OFS='\t' \
    -v fai="$REF_FAI" -v genes="$GENE_LIST" \
    -v out="$OUT_TSV" -v normout="$NORM_TSV" '
  BEGIN {
    # 参考 FASTA 的 .fai：第 1 列转录本 ID，第 2 列序列长度（bp）
    while ((getline line < fai) > 0) {
      split(line, a, "\t")
      len[a[1]] = a[2] + 0
    }
    close(fai)
    # 保存基因顺序，保证输出可复现
    while ((getline line < genes) > 0) {
      keep[line] = 1
      gene_order[++ngene] = line
    }
    close(genes)
    print "genotype", "sample", "feature_id", "transcript_length_bp",
          "tpm_denominator_count_per_kb", "raw_count", "tpm" > out
    print "sample", "genotype", "tpm_denominator_count_per_kb" > normout
  }
  # count 矩阵表头：只保留目标基因型 27°C 的 rep3/4/5
  FNR == 1 {
    for (i = 2; i <= NF; i++) {
      sample = $i
      if (sample ~ /^Mock-27-(WT|trua1b_1|trua1b_2)\.rep[345]$/) {
        ns++
        col[ns] = i
        samples[ns] = sample
        geno = sample
        sub(/^Mock-27-/, "", geno)
        sub(/\.rep[345]$/, "", geno)
        genotypes[ns] = geno
      }
    }
    if (ns != 9) {
      printf "[data_prep] ERROR: expected 9 target samples, found %d\n", ns > "/dev/stderr"
      exit 1
    }
    next
  }
  {
    if (!($1 in len) || len[$1] <= 0) {
      printf "[data_prep] ERROR: missing/invalid transcript length for %s\n", $1 > "/dev/stderr"
      exit 1
    }
    # TPM 分母必须累加全部转录本的长度校正丰度，不能只累加上调基因
    for (j = 1; j <= ns; j++) {
      value = $col[j] + 0
      rate = value / (len[$1] / 1000)
      denom[j] += rate
      if ($1 in keep) raw[$1, j] = value
    }
  }
  END {
    for (j = 1; j <= ns; j++) {
      if (denom[j] <= 0) {
        printf "[data_prep] ERROR: nonpositive TPM denominator for %s\n", samples[j] > "/dev/stderr"
        exit 1
      }
      print samples[j], genotypes[j], sprintf("%.6f", denom[j]) > normout
    }
    for (j = 1; j <= ns; j++) {
      for (i = 1; i <= ngene; i++) {
        id = gene_order[i]
        tpm = (raw[id, j] / (len[id] / 1000)) / denom[j] * 1000000
        print genotypes[j], samples[j], id, len[id], sprintf("%.6f", denom[j]),
              sprintf("%.0f", raw[id, j]), sprintf("%.6f", tpm) > out
      }
    }
  }
' "$COUNT_TSV"

N_ROWS=$(( $(wc -l < "$OUT_TSV") - 1 ))
EXPECTED_ROWS=$(( N_GENES * 9 ))
[ "$N_ROWS" -eq "$EXPECTED_ROWS" ] || {
  echo "[data_prep] ERROR: expected $EXPECTED_ROWS rows, got $N_ROWS" >&2
  exit 1
}
echo "[data_prep] WT 27C-up genes: $N_GENES; TPM rows: $N_ROWS -> $OUT_TSV"
echo "[data_prep] TPM normalization denominators: $NORM_TSV"
