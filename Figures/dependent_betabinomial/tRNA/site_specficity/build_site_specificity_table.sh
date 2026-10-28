#!/usr/bin/env bash
# build_site_specificity_table.sh — 位点特异性（site specificity）分析的数据构建
#
# 问题：WT 里在某个 canonical 位号（Sprinzl/fixed 后）上有修饰的代表 tRNA 有多少个？
#       这些 tRNA 在一组目标突变体（各两条株系）里修饰比例变成什么样？
#
# 入选口径（四条件同时满足）：
#   1. 位点在 WT 的 Fisher positive 宇宙里（WT-<tissue>_sites.Treat_avg_ratio.tsv，
#      表内全部行 Pass=1，即 WT 三个 Treat 重复都检出 CMC 效应 = WT 中有修饰）；
#   2. Sprinzl_pos == 指定位号（native pos 各异，26/27/28/29 都可能映射到 27）；
#   3. chr 是 family_top 的代表 tRNA（每个 family mean_CPM 最高的成员）；
#   4. 分室前缀在 compartments 白名单内（空 = 全部分室）。
#
# 比例口径：WT 取 Treat_avg_ratio；每条株系取 <g>_<line>-<tissue>.dependent-sites.tsv
#   的 Mut_Ratio（该株系 Treat 重复的均值，未筛选），与基因型图组一致；
#   株系表缺该位点时留 NA（绝不当作 0% 修饰），同时带出该株系的 per-line is_dependent。
#   WT 表在各基因型目录下逐字节相同，取第一个基因型目录的即可。
#
# 同一 ncRNA 参考混有 sn/snoRNA，按 ^tRNA- 筛掉（与 build_genotype_dependent_tables.sh 同口径）。
#
# Usage: bash build_site_specificity_table.sh <tissue> <genotypes> <sprinzl_pos> <output_dir> [compartments]
#   tissue ∈ Root|leaf|Flower；genotypes 逗号分隔（如 trua1a,trua5）；
#   sprinzl_pos 如 27；compartments 逗号分隔的 chr 前缀分室（Nu/Mt/Pt），留空 = 全部。
set -euo pipefail
IFS=$'\n\t'

if [[ $# -lt 4 ]]; then
  echo "Usage: $0 <tissue> <genotypes> <sprinzl_pos> <output_dir> [compartments]" >&2
  exit 1
fi

readonly TISSUE="$1"
readonly GENOTYPES="$2"
readonly SPRINZL_POS="$3"
readonly OUT_DIR="$4"
readonly COMPARTMENTS="${5:-}"

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# 本目录比基因型目录（build_genotype_dependent_tables.sh 所在层）深一级，故多退一层
readonly ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
readonly DEPENDENT_DIR="$ROOT_DIR/tRNA/STAT/dependent_betabinomial/ncRNA"
readonly REP_FILE="$ROOT_DIR/tRNA/count/ncRNA/family_top/family_top_tRNA.tsv"
readonly SPRINZL_SCRIPT="$ROOT_DIR/BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py"

readonly WT_FILE="$DEPENDENT_DIR/$TISSUE/${GENOTYPES%%,*}/WT-${TISSUE}_sites.Treat_avg_ratio.tsv"

log() { printf '[site_specificity] %s\n' "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

[[ -s "$WT_FILE" ]] || error "Missing WT universe: $WT_FILE"
[[ -s "$REP_FILE" ]] || error "Missing representative table: $REP_FILE"
[[ -s "$SPRINZL_SCRIPT" ]] || error "Missing annotation script: $SPRINZL_SCRIPT"
mkdir -p "$OUT_DIR"

IFS=',' read -r -a GENOTYPE_LIST <<< "$GENOTYPES"
[[ ${#GENOTYPE_LIST[@]} -ge 1 ]] || error "genotypes must be a non-empty comma list"

# 每基因型两份株系表；同时核对各基因型目录下的 WT 表与取用的那份逐字节一致
LINE_FILES=()
for g in "${GENOTYPE_LIST[@]}"; do
  for f in "$DEPENDENT_DIR/$TISSUE/$g/${g}_1-${TISSUE}.dependent-sites.tsv" \
           "$DEPENDENT_DIR/$TISSUE/$g/${g}_2-${TISSUE}.dependent-sites.tsv"; do
    [[ -s "$f" ]] || error "Missing line-level table: $f"
    LINE_FILES+=("$f")
  done
  wt_other="$DEPENDENT_DIR/$TISSUE/$g/WT-${TISSUE}_sites.Treat_avg_ratio.tsv"
  if [[ -s "$wt_other" ]] && ! cmp -s "$WT_FILE" "$wt_other"; then
    error "WT universe differs across genotype dirs: $WT_FILE vs $wt_other"
  fi
done

# WT positive 宇宙原表不带 Sprinzl_pos（注释只对 21 张 final dependent 表做过），
# 在共享层（本脚本所在目录）补一份注释版中间产物：同组织的多个位号文件夹共用，
# 不动 tRNA/STAT/ 原始数据。
readonly WT_SPRINZL_FILE="$SCRIPT_DIR/WT-${TISSUE}_sites.Sprinzl.tsv"
readonly STEM_TISSUE="$(printf '%s' "$TISSUE" | tr '[:upper:]' '[:lower:]')"
readonly OUT_TSV="$OUT_DIR/pos${SPRINZL_POS}_${STEM_TISSUE}_site_table.tsv"

python3 "$SPRINZL_SCRIPT" -i "$WT_FILE" -o "$WT_SPRINZL_FILE" >&2

# 表头：株系值两列 + per-line dependent 两列，按基因型顺序展开
{
  printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tTreat_avg_depth'
  for g in "${GENOTYPE_LIST[@]}"; do
    printf '\t%s_1\t%s_2\t%s_1_dependent\t%s_2_dependent' "$g" "$g" "$g" "$g"
  done
  printf '\n'
} > "$OUT_TSV"

# 株系表路径列表传给 awk，按 FILENAME 在列表中的序号分派到 (基因型, 株系)
LINE_FILE_ARG="$(IFS=','; printf '%s' "${LINE_FILES[*]}")"

awk -F'\t' -v OFS='\t' \
    -v target_pos="$SPRINZL_POS" \
    -v genotypes="$GENOTYPES" \
    -v tissue="$TISSUE" \
    -v compartments="$COMPARTMENTS" \
    -v trna_pat='^tRNA-' \
    -v rep_file="$REP_FILE" \
    -v wt_file="$WT_SPRINZL_FILE" \
    -v line_file_arg="$LINE_FILE_ARG" '
# 列名 -> 列号（各输入列序不保证一致，一律按表头取）
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
function num_ok(x) { return x ~ /^-?[0-9]+(\.[0-9]+)?([eE][-+]?[0-9]+)?$/ }
FILENAME == rep_file {
    if (FNR == 1) { load_cols(); require_cols("feature_id"); next }
    rep[$(col["feature_id"])] = 1
    next
}
FILENAME == wt_file {
    if (FNR == 1) { load_cols(); require_cols("chr,pos,Sprinzl_pos,strand,context,Treat_avg_ratio,Treat_avg_depth,Pass"); next }
    if ($1 !~ trna_pat) next
    if (!($1 in rep)) next
    # 分室限定（空 = 不过滤）：chr 前缀 tRNA-Nu-/tRNA-Mt-/tRNA-Pt- 对应 Nu/Mt/Pt
    if (compartments != "") {
        comp_ok = 0
        n_comp = split(compartments, comp_list, ",")
        for (ci = 1; ci <= n_comp; ci++) {
            if (index($1, "tRNA-" comp_list[ci] "-") == 1) { comp_ok = 1; break }
        }
        if (!comp_ok) next
    }
    if ($(col["Sprinzl_pos"]) != target_pos) next
    if ($(col["Pass"]) != "1") {
        printf("WT site not Pass=1 in positive universe: %s %s\n", $1, $2) > "/dev/stderr"
        exit 1
    }
    if (!num_ok($(col["Treat_avg_ratio"]))) {
        printf("Non-numeric Treat_avg_ratio for %s %s\n", $1, $2) > "/dev/stderr"
        exit 1
    }
    skey = $1 SUBSEP $2 SUBSEP $(col["strand"])
    if (skey in wt) {
        printf("Duplicate WT site after filtering: %s %s\n", $1, $2) > "/dev/stderr"
        exit 1
    }
    wt[skey] = $(col["Treat_avg_ratio"]) + 0
    wt_depth[skey] = $(col["Treat_avg_depth"])
    ctx[skey] = $(col["context"])
    order[++n] = skey
    next
}
# 株系表：FILENAME 在 line_file_arg 里的序号 -> (基因型序号 gi, 株系 ln)
FNR == 1 {
    load_cols()
    require_cols("chr,pos,strand,context,WT_Ratio,Mut_Ratio,is_dependent")
    n_files = split(line_file_arg, lf, ",")
    file_rank = 0
    for (fi = 1; fi <= n_files; fi++) {
        if (lf[fi] == FILENAME) { file_rank = fi; break }
    }
    if (file_rank == 0) {
        printf("Unexpected input file (not in line list): %s\n", FILENAME) > "/dev/stderr"
        exit 1
    }
    split(genotypes, glist, ",")
    gi = int((file_rank - 1) / 2) + 1
    ln = ((file_rank - 1) % 2) + 1
    next
}
{
    skey = $1 SUBSEP $2 SUBSEP $(col["strand"])
    if (!(skey in wt)) next
    wt_now = $(col["WT_Ratio"]) + 0
    if ((wt[skey] - wt_now > 1e-6) || (wt_now - wt[skey] > 1e-6)) {
        printf("WT ratio mismatch between WT universe and %s: %s %s\n", FILENAME, $1, $2) > "/dev/stderr"
        exit 1
    }
    if (ctx[skey] != $(col["context"])) {
        printf("Context mismatch between WT universe and %s: %s %s\n", FILENAME, $1, $2) > "/dev/stderr"
        exit 1
    }
    mut[gi, ln, skey] = $(col["Mut_Ratio"])
    dep[gi, ln, skey] = $(col["is_dependent"])
}
END {
    split(genotypes, glist, ",")
    n_g = length(glist)
    for (i = 1; i <= n; i++) {
        skey = order[i]
        split(skey, f, SUBSEP)
        row = f[1] OFS f[2] OFS target_pos OFS f[3] OFS ctx[skey] OFS \
              sprintf("%.15g", wt[skey]) OFS wt_depth[skey]
        for (gi = 1; gi <= n_g; gi++) {
            # 列序与表头一致：两株系值在前，两株系 per-line dependent 在后
            for (ln = 1; ln <= 2; ln++) {
                row = row OFS ((gi, ln, skey) in mut ? mut[gi, ln, skey] : "NA")
            }
            for (ln = 1; ln <= 2; ln++) {
                row = row OFS ((gi, ln, skey) in mut ? dep[gi, ln, skey] : "NA")
            }
        }
        print row
    }
    for (gi = 1; gi <= n_g; gi++) {
        miss[gi] = 0
        for (i = 1; i <= n; i++) {
            for (ln = 1; ln <= 2; ln++) {
                if (!((gi, ln, order[i]) in mut)) miss[gi]++
            }
        }
    }
    msg = ""
    for (gi = 1; gi <= n_g; gi++) {
        msg = msg sprintf("%s%s: %d NA", (gi > 1 ? " / " : ""), glist[gi], miss[gi])
    }
    printf("[site_specificity] %s pos%s%s: %d selected tRNA sites (%s)\n", \
           tissue, target_pos, (compartments == "" ? " all compartments" : " compartments=" compartments), \
           n, msg) > "/dev/stderr"
    if (n == 0) {
        printf("No representative tRNA sites selected for pos%s.\n", target_pos) > "/dev/stderr"
        exit 1
    }
}
' "$REP_FILE" "$WT_SPRINZL_FILE" "${LINE_FILES[@]}" >> "$OUT_TSV"

log "Generated: $OUT_TSV ($(( $(wc -l < "$OUT_TSV") - 1 )) sites, genotypes: $GENOTYPES)"
