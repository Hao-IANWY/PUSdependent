#!/usr/bin/env bash

# 构建每个 WT 组织的"表达 U 位点池"(expressed pool)并按 mRNA 区域注释。
#
# 表达池 = 该组织样本中实际被测序覆盖到、具备修饰检测能力的尿苷位点全集, 作为区域修饰率
# 的分母。背景: polyARNA/STAT/fisher/mRNA/WT-*_sites.Fisher.tsv 并非全部进入检验的位点 ——
# Pseudo_callsites.R 在写出前已过滤 pooled Fisher p.adj < 0.01。表达池按同一脚本的
# 入池规则从每重复候选文件重建:
#   1. STAT/adjustGap/mRNA/*.filter_T.tsv: 每重复中参考碱基为 T 且 depth >= 10 的位点;
#   2. Mock x3 + Treat x3 共 6 个重复按位点取并集, 缺失重复的 depth 记 0;
#   3. Mock 组 depth 总和 >= 30 且 Treat 组总和 >= 30 (即组平均 >= min_avg_depth = 10)。
# 池位点随后走与阳性位点完全相同的注释流程(gppy t2g 转基因组坐标 + bed_annoat.R 区域注释),
# 保证分子(阳性位点)与分母(表达池)的区域口径一致。
#
# 输出: site_analyse/WT-<Tissue>.expressed_pool.annoated.tsv (chr/pos/strand/annoation)

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

    # 并集 + 组均 depth >= 10 (Mock/Treat 各 3 重复, 缺失记 0, 总和阈值 = 3 x 10)
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

    # 行号作为 id 追加为第 4 列(t2g 只读前两列 chr/pos): t2g 映射失败(转录本不在 GTF)的行
    # 在 pool 与 bed 两侧按 id 一致剔除, 保证与注释输出的行对齐
    local numbered_tsv="${OUTPUT_DIR}/${prefix}.expressed_pool.numbered.tsv"
    awk 'BEGIN{OFS="\t"}{print $0, NR}' "$pool_tsv" > "$numbered_tsv"

    # t2g 对百万级行数较慢, 分块并行; 各块输出按块序拼接以保持与 pool.tsv 行序一致
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

    # 中间文件(pool/mapped/bed/annoated.bed)可由本脚本完全再生, 只保留最终注释表
    rm -f "$pool_tsv" "$mapped_tsv" "$pool_bed" "$annotated_bed"

    log "Wrote ${annotated_tsv}"
  done

  log "All tested pools built and annotated. Outputs are in ${OUTPUT_DIR}"
}

main "$@"
