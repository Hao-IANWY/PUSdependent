#!/usr/bin/env bash
# build_genotype_dependent_tables.sh — 汇总单个 trua 突变体的 tRNA dependent 位点表
#
# 移植自 archive/20260902/paper/tRNA_dependent/build_union_dependent_sites.sh，
# 数据源从旧 Fisher 口径换成 tRNA/STAT/dependent_methylKit（methylKit 口径）。
#
# 7 个基因型的逻辑完全一致，只有 GENOTYPE 不同，故收在一个共享脚本里，
# 各 <genotype>/data_prep.sh 只做三行转发。
#
# 与旧口径的三处差异：
#   1. dependent 表（<g>-<tissue>.dependent.tsv）没有 fix_pos。canonical 位号在这里现算：
#      按 All-tRNA_rep-structure.csv 的区段顺序/固定长度把 native pos 重标定到 0-76。
#   2. dependent 表也没有每株系修饰比。从 <g>_{1,2}-<tissue>.dependent-sites.tsv
#      的 Mut_Ratio（该株系 3 个 Treat 重复的均值，未筛选）按 chr+pos+strand 取。
#      dependent 表本身就是两株系 AND 的结果，所以必定取得到。
#   3. 同一 ncRNA 参考里混有 sn/snoRNA（589 条 = 220 tRNA + 287 snoRNA + 82 snRNA），
#      dependent 表里也有它们的位点，必须按 chr 前缀 ^tRNA- 筛掉。
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
readonly DEPENDENT_DIR="$ROOT_DIR/tRNA/STAT/dependent_methylKit/ncRNA"
readonly REP_STRUCTURE_CSV="$SCRIPT_DIR/All-tRNA_rep-structure.csv"

readonly LONG_OUT="$OUT_DIR/${GENOTYPE}_all_tissues.dependent.tsv"
readonly UNION_OUT="$OUT_DIR/${GENOTYPE}_union.dependent.tsv"
readonly TRNA_CHR_PATTERN='^tRNA-'

readonly -a TISSUES=("Root" "leaf" "Flower")
readonly -a LINES=("1" "2")

log() { printf '[data_prep] %s\n' "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

[[ -s "$REP_STRUCTURE_CSV" ]] || error "Missing representative structure CSV: $REP_STRUCTURE_CSV"
[[ -d "$DEPENDENT_DIR" ]] || error "Missing dependent directory: $DEPENDENT_DIR"
mkdir -p "$OUT_DIR"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------------------------------------------------------------------------
# 1. native pos -> canonical 位号（与 archive dependent_site_overview/data.sh 一致）
#    区段顺序 / 固定长度 / 重标定公式逐字沿用；round 半值取上（archive 的 awk 行为）。
#    输出 <rep_id>\t<native_pos>\t<canonical_pos>
# ---------------------------------------------------------------------------
CANONICAL_MAP="$TMP_DIR/canonical-pos-map.tsv"
awk -F',' 'BEGIN {
    OFS = "\t"
    n_regions = split("p1,Accstem_5p,p8-9,Dstem_5p,Dloop,Dstem_3p,p26,Acstem_5p,Acloop,Acstem_3p,Vregion,Tstem_5p,Tloop,Tstem_3p,Accstem_3p,p73,CCA", region_order, ",")
    split("1,7,2,4,8,4,1,5,7,5,5,5,7,5,7,1,3", fixed_lengths, ",")

    canonical_offset = 0
    for (i = 1; i <= n_regions; i++) {
      region = region_order[i]
      fixed_len_map[region] = fixed_lengths[i] + 0
      region_offset[region] = canonical_offset
      canonical_offset += fixed_len_map[region]
    }
  }
  NR == 1 {
    for (i = 3; i <= NF; i++) {
      gsub(/\r/, "", $i)
      region_name[i] = $i
    }
    next
  }
  {
    rep_id = $2
    gsub(/\r/, "", rep_id)

    native_start = 1
    for (i = 3; i <= NF; i++) {
      seq = $i
      gsub(/\r/, "", seq)

      region = region_name[i]
      region_len = length(seq)
      if (region_len == 0) {
        continue
      }

      current_fixed_len = fixed_len_map[region]
      if (current_fixed_len == 0) {
        printf "Unknown structural region: %s\n", region > "/dev/stderr"
        exit 1
      }

      for (native_pos = 1; native_pos <= region_len; native_pos++) {
        if (current_fixed_len == 1 || region_len == 1) {
          canonical_region_pos = 1
        } else {
          canonical_region_pos = int(((native_pos - 1) * (current_fixed_len - 1) / (region_len - 1)) + 0.5) + 1
        }

        canonical_pos = region_offset[region] + canonical_region_pos - 1
        print rep_id, native_start + native_pos - 1, canonical_pos
      }

      native_start += region_len
    }
  }' "$REP_STRUCTURE_CSV" > "$CANONICAL_MAP"

test -s "$CANONICAL_MAP" || error "Failed to build canonical position map from $REP_STRUCTURE_CSV"

# ---------------------------------------------------------------------------
# 2. 三组织纵向拼接：dependent 表 + 每株系 Mut_Ratio + fix_pos + Tissue
# ---------------------------------------------------------------------------
TMP_LONG="$TMP_DIR/all-tissues.tsv"
printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tMut_Ratio\tDelta_Ratio(WT-Mut)\t%s\t%s\tTissue\n' \
  "${GENOTYPE}_1" "${GENOTYPE}_2" > "$TMP_LONG"

for tissue in "${TISSUES[@]}"; do
  final_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/$GENOTYPE-$tissue.dependent.tsv"
  line1_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_1-$tissue.dependent-sites.tsv"
  line2_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_2-$tissue.dependent-sites.tsv"

  [[ -s "$final_file" ]] || error "Missing dependent table: $final_file"
  [[ -s "$line1_file" ]] || error "Missing line-level table: $line1_file"
  [[ -s "$line2_file" ]] || error "Missing line-level table: $line2_file"

  awk -F'\t' -v OFS='\t' \
      -v tissue="$tissue" \
      -v trna_pat="$TRNA_CHR_PATTERN" \
      -v map_file="$CANONICAL_MAP" \
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
  FILENAME == map_file {
      canonical[$1 SUBSEP $2] = $3
      next
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
  # 第三个位置上的 dependent 表
  FNR == 1 { load_cols(); require_cols("chr,pos,strand,context,WT_Ratio,Mut_Ratio,Delta_Ratio(WT-Mut)"); next }
  {
      chr = $(col["chr"])
      # 只保留 tRNA 来源位点，丢掉同一参考里的 sn/snoRNA
      if (chr !~ trna_pat) next

      pos = $(col["pos"])
      strand = $(col["strand"])
      rep_id = chr
      sub(/^tRNA-/, "", rep_id)
      sub(/-[0-9]+$/, "", rep_id)

      mkey = rep_id SUBSEP pos
      if (!(mkey in canonical)) {
          printf "Failed to map %s native position %s in %s\n", rep_id, pos, FILENAME > "/dev/stderr"
          exit 1
      }

      skey = chr SUBSEP pos SUBSEP strand
      print chr, pos, canonical[mkey], strand, $(col["context"]), \
            $(col["WT_Ratio"]), $(col["Mut_Ratio"]), $(col["Delta_Ratio(WT-Mut)"]), \
            (skey in mut1 ? mut1[skey] : 0), (skey in mut2 ? mut2[skey] : 0), tissue
  }' "$CANONICAL_MAP" "$line1_file" "$line2_file" "$final_file" >> "$TMP_LONG"

  log "$GENOTYPE / $tissue done"
done

cp "$TMP_LONG" "$LONG_OUT"

# ---------------------------------------------------------------------------
# 3. 跨组织并集：每个位点取 Delta 最大的那次观测，附 Tissues / n_tissues / peak_tissue
#    沿用 archive build_union_dependent_sites.sh 的键与排序。
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
