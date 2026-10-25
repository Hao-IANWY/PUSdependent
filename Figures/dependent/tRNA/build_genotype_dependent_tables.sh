#!/usr/bin/env bash
# build_genotype_dependent_tables.sh — 汇总单个 trua 突变体的 tRNA dependent 位点表
#
# 数据源：tRNA/STAT/dependent/ncRNA（纯效应量 Delta>=0.05 口径，无显著性检验），且已由
#   BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py
#   加上 Sprinzl 位号（文件名 <g>-<tissue>.dependent.Sprinzl.tsv）。
#
# 与旧版的三处差异：
#   1. fix_pos 直接取上游的 Sprinzl_pos —— 本目录不再做任何坐标转换，
#      也不再依赖 All-tRNA_rep-structure.csv。
#   2. 上游按每个拷贝自己的区段注释计算，长度/区段划分与家族代表不同的拷贝也正确。
#   3. 每株系修饰比仍从 <g>_{1,2}-<tissue>.dependent-sites.tsv 的 Mut_Ratio
#      （该株系 3 个 Treat 重复的均值，未筛选）按 chr+pos+strand 取；这些文件不需要注释列。
#
# 同一 ncRNA 参考里混有 sn/snoRNA，按 chr 前缀 ^tRNA- 筛掉。
#
# Usage: bash build_genotype_dependent_tables.sh <genotype> <output_dir>
set -euo pipefail
IFS=$'\n\t'

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <genotype> <output_dir>" >&2
  exit 1
fi

readonly GENOTYPE="$1"
readonly OUT_DIR="$2"

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "$SCRIPT_DIR/../../.." && pwd)"
readonly DEPENDENT_DIR="$ROOT_DIR/tRNA/STAT/dependent/ncRNA"

readonly LONG_OUT="$OUT_DIR/${GENOTYPE}_all_tissues.dependent.tsv"
readonly UNION_OUT="$OUT_DIR/${GENOTYPE}_union.dependent.tsv"
readonly TRNA_CHR_PATTERN='^tRNA-'

readonly -a TISSUES=("Root" "leaf" "Flower")
readonly -a LINES=("1" "2")

log() { printf '[data_prep] %s\n' "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

[[ -d "$DEPENDENT_DIR" ]] || error "Missing dependent directory: $DEPENDENT_DIR"
mkdir -p "$OUT_DIR"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------------------------------------------------------------------------
# 三组织纵向拼接：dependent 表（自带 Sprinzl_pos）+ 每株系 Mut_Ratio + Tissue
# ---------------------------------------------------------------------------
TMP_LONG="$TMP_DIR/all-tissues.tsv"
printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tMut_Ratio\tDelta_Ratio(WT-Mut)\t%s\t%s\tTissue\n' \
  "${GENOTYPE}_1" "${GENOTYPE}_2" > "$TMP_LONG"

for tissue in "${TISSUES[@]}"; do
  final_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/$GENOTYPE-$tissue.dependent.Sprinzl.tsv"
  line1_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_1-$tissue.dependent-sites.tsv"
  line2_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_2-$tissue.dependent-sites.tsv"

  [[ -s "$final_file" ]] || error "Missing annotated dependent table: $final_file（先跑 add_sprinzl_pos.py）"
  [[ -s "$line1_file" ]] || error "Missing line-level table: $line1_file"
  [[ -s "$line2_file" ]] || error "Missing line-level table: $line2_file"

  awk -F'\t' -v OFS='\t' \
      -v tissue="$tissue" \
      -v trna_pat="$TRNA_CHR_PATTERN" \
      -v line1_file="$line1_file" \
      -v line2_file="$line2_file" '
  # 列名 -> 列号（三份输入的列序不保证一致，一律按表头取）
  function load_cols(    i, name) {
      delete col
      for (i = 1; i <= NF; i++) {
          name = $i
          sub(/\r$/, "", name)
          col[name] = i
      }
  }
  function require_cols(required,    i, fields, name) {
      split(required, fields, ",")
      for (i in fields) {
          name = fields[i]
          if (!(name in col)) {
              printf("Required column missing in %s: %s\n", FILENAME, name) > "/dev/stderr"
              exit 1
          }
      }
  }
  FILENAME == line1_file {
      if (FNR == 1) { load_cols(); require_cols("chr,pos,strand,Mut_Ratio"); next }
      mut1[$1 SUBSEP $2 SUBSEP $3] = $(col["Mut_Ratio"])
      next
  }
  FILENAME == line2_file {
      if (FNR == 1) { load_cols(); require_cols("chr,pos,strand,Mut_Ratio"); next }
      mut2[$1 SUBSEP $2 SUBSEP $3] = $(col["Mut_Ratio"])
      next
  }
  # 第三个位置上的 dependent 表（带 Sprinzl_pos）
  FNR == 1 { load_cols(); require_cols("chr,pos,Sprinzl_pos,strand,context,WT_Ratio,Mut_Ratio,Delta_Ratio(WT-Mut)"); next }
  {
      chr = $(col["chr"])
      if (chr !~ trna_pat) next                  # 丢掉同一参考里的 sn/snoRNA

      pos = $(col["pos"])
      strand = $(col["strand"])
      fix_pos = $(col["Sprinzl_pos"])
      if (fix_pos == "" || fix_pos !~ /^[0-9]+$/ || fix_pos + 0 < 0 || fix_pos + 0 > 76) {
          printf("Missing/invalid Sprinzl_pos for %s pos %s in %s\n", chr, pos, FILENAME) > "/dev/stderr"
          exit 1
      }

      skey = chr SUBSEP pos SUBSEP strand
      print chr, pos, fix_pos + 0, strand, $(col["context"]), \
            $(col["WT_Ratio"]), $(col["Mut_Ratio"]), $(col["Delta_Ratio(WT-Mut)"]), \
            (skey in mut1 ? mut1[skey] : 0), (skey in mut2 ? mut2[skey] : 0), tissue
  }' "$line1_file" "$line2_file" "$final_file" >> "$TMP_LONG"

  log "$GENOTYPE / $tissue done"
done

cp "$TMP_LONG" "$LONG_OUT"

# ---------------------------------------------------------------------------
# 跨组织并集：每个位点取 Delta 最大的那次观测，附 Tissues / n_tissues / peak_tissue
# ---------------------------------------------------------------------------
TMP_BODY="$TMP_DIR/union-body.tsv"
awk 'BEGIN {
    FS = OFS = "\t"
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) {
      idx[$i] = i
    }
    next
  }
  {
    key = $1 OFS $2 OFS $3 OFS $4 OFS $5
    tissue = $idx["Tissue"]
    wt_ratio = $idx["WT_Ratio"] + 0
    mut_ratio = $idx["Mut_Ratio"] + 0
    delta_ratio = $idx["Delta_Ratio(WT-Mut)"] + 0

    if (!(key in best_delta) || delta_ratio > best_delta[key]) {
      best_delta[key] = delta_ratio
      best_wt[key] = wt_ratio
      best_mut[key] = mut_ratio
      peak_tissue[key] = tissue
    }

    if (!(seen_tissue[key SUBSEP tissue]++)) {
      tissues[key] = tissues[key] ? tissues[key] "," tissue : tissue
      tissue_count[key]++
    }

    row[key] = $1 OFS $2 OFS $3 OFS $4 OFS $5
  }
  END {
    for (key in row) {
      print row[key], best_wt[key], best_mut[key], best_delta[key], tissues[key], tissue_count[key], peak_tissue[key]
    }
  }' "$TMP_LONG" > "$TMP_BODY"

{
  printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tMut_Ratio\tDelta_Ratio(WT-Mut)\tTissues\tn_tissues\tpeak_tissue\n'
  sort -t $'\t' -k1,1 -k3,3n -k2,2n "$TMP_BODY"
} > "$UNION_OUT"

log "$GENOTYPE: $(( $(wc -l < "$LONG_OUT") - 1 )) tissue-resolved rows, $(( $(wc -l < "$UNION_OUT") - 1 )) union sites"
