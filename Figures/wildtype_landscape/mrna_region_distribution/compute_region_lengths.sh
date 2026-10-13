#!/usr/bin/env bash

# 计算 mRNA 参考(代表基因模型)各区域的拼接序列总长, 作为 "每 1000 nt 修饰密度" 图的分母。
#
# mRNA 参考为 Araport11 representative gene model(每基因 1 个转录本), 无多 isoform 重复计数
# 问题; 区域长度来自 gppy txinfo(基于外显子结构的拼接长度, utr3 = tx_len - cds_len - utr5_len),
# 并用参考 .fai 的 length 列总和交叉校验 sum(tx_len)。
#
# 输出: mrna_region_distribution_region_lengths.tsv (Region/Length_nt, 5'UTR/CDS/3'UTR 三行)

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

  # 不用 local: EXIT trap 在 main 返回后才触发, local 变量已出作用域
  tmp_dir="$(mktemp -d "${SCRIPT_DIR}/.region_lengths.XXXXXX")"
  trap 'rm -rf "$tmp_dir"' EXIT

  # 参考转录本 id 集 + 拼接总长(.fai 第 2 列)
  cut -f1 "${MRNA_REF}.fai" > "$tmp_dir/ref_ids.txt"
  local n_ref fai_total
  n_ref="$(wc -l < "$tmp_dir/ref_ids.txt")"
  fai_total="$(awk '{s += $2} END {print s}' "${MRNA_REF}.fai")"
  log "Reference transcripts: ${n_ref}, total spliced length (.fai): ${fai_total}"

  # 每个转录本的拼接区域长度: tx_name gene_id chrom strand nexon tx_len cds_len utr5_len utr3_len
  # (用 txinfo_basic: txinfo 解析属性时对细胞器基因报 KeyError)
  "$PYTHON_BIN" "$GPPY" txinfo_basic -g "$GTF" > "$tmp_dir/tx_info.tsv"

  # 过滤到参考集合并按区域求和; utr3_len 为 0 的非编码转录本自然不贡献 CDS
  awk -F'\t' '
    NR==FNR { ref[$1]=1; next }
    FNR==1 { next }
    ($1 in ref) { tx+=$6; cds+=$7; u5+=$8; u3+=$9; n++ }
    END {
      # 校验行(matched_tx/sum_tx_len)走 stderr 落盘, 区域长度行走 stdout
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

  # 列出参考中未匹配 GTF 的 id(如细胞器转录本), 供核查
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
