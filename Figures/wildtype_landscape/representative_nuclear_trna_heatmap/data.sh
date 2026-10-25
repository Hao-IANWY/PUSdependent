#!/usr/bin/env bash
# data.sh — 把 WT 各组织的 tRNA positive site 折叠成
#           代表性核 tRNA × canonical 0–76 位 的 ratio 表，供 plot.R 画热图。
#
# 与旧版的三处区别：
#   1. 代表 tRNA 取自 tRNA/count/ncRNA/family_top/family_top_tRNA.tsv
#      （每个 anticodon 家族里 mean_CPM 最高的成员，即表达量口径的代表）。
#   2. Sprinzl 坐标不再在本目录转换：直接读 tRNA/STAT/fisher/ncRNA 下已注释好的
#      WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv 的 Sprinzl_pos 列。
#      每个拷贝用自己的区段注释，因此长度与家族代表不一致的拷贝也不再剔除。
#   3. 行仍是家族层（同 family 的所有拷贝位点都并入同一行）；
#      代表只决定该行的 Representative_ID 标签。
#
# 输入：
#   tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.Sprinzl.tsv
#       chr, pos, Sprinzl_pos, strand, ratio, depth
#   tRNA/count/ncRNA/family_top/family_top_tRNA.tsv
#       feature_id, family, family_size, mean_CPM, ...
#
# 输出：
#   representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv
#       tissue, Rep_ID, Representative_ID, chr, pos, strand, ratio, depth, canonical_pos
#   trna_landscape_summary.tsv
#       tissue, sites_in, excluded_nonNu, unmapped, sites_mapped, reps_covered, cells
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
readonly INPUT_DIR="$REPO_ROOT/tRNA/STAT/fisher/ncRNA"
readonly REPS="$REPO_ROOT/tRNA/count/ncRNA/family_top/family_top_tRNA.tsv"
readonly COLUMNS=$'tissue\tRep_ID\tRepresentative_ID\tchr\tpos\tstrand\tratio\tdepth\tcanonical_pos'

log() { printf '[data_prep] %s\n' "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

[[ -s "$REPS" ]] || error "Missing family_top table: $REPS"

readonly SUMMARY="$SCRIPT_DIR/trna_landscape_summary.tsv"
readonly COUNTS="$SCRIPT_DIR/.counts.tmp"
trap 'rm -f "$COUNTS"' EXIT

printf 'tissue\tsites_in\texcluded_nonNu\tunmapped\tsites_mapped\treps_covered\tcells\n' > "$SUMMARY"

for tissue in Root Leaf Flower; do
    input="$INPUT_DIR/WT_${tissue}_sites.positive.Sprinzl.tsv"
    [[ -s "$input" ]] || error "Missing annotated input: $input（先跑 add_sprinzl_pos.py）"
    prefix="$(printf '%s' "$tissue" | tr '[:upper:]' '[:lower:]')"
    out="$SCRIPT_DIR/representative_nuclear_trna_mapped_sites_${prefix}.tsv"

    printf '%s\n' "$COLUMNS" > "$out"
    awk -F'\t' -v OFS='\t' -v tissue="$tissue" -v reps="$REPS" -v counts="$COUNTS" '
    # 家族代表表：只取核 tRNA，family 名即 Rep_ID（如 Nu-Arg-TCG）
    FILENAME == reps {
        if (FNR > 1 && $1 ~ /^tRNA-Nu-/) REP[$2] = $1
        next
    }
    FNR == 1 { next }
    {
        sites_in++
        chr = $1
        if (chr !~ /^tRNA-Nu-/) { excluded++; next }        # 只画核 tRNA
        fam = chr
        sub(/^tRNA-Nu-/, "Nu-", fam)
        sub(/-[0-9]+$/, "", fam)
        if (!(fam in REP)) { excluded++; next }             # family_top 里没有的家族

        canon = $3
        if (canon == "" || canon !~ /^[0-9]+$/ || canon + 0 < 0 || canon + 0 > 76) {
            unmapped++; next                                # Sprinzl_pos 缺失/越界
        }

        print tissue, fam, REP[fam], chr, $2, $4, sprintf("%.10g", $5), sprintf("%.10g", $6), canon + 0
        mapped++
        seen_rep[fam] = 1
        seen_cell[fam SUBSEP canon] = 1
    }
    END {
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", tissue, sites_in + 0, excluded + 0,
               unmapped + 0, mapped + 0, length(seen_rep), length(seen_cell) > counts
    }
    ' "$REPS" "$input" | sort -t$'\t' -k2,2 -k9,9n -k5,5n >> "$out"

    counts_line="$(cat "$COUNTS")"
    printf '%s\n' "$counts_line" >> "$SUMMARY"
    log "$(printf '%s' "$counts_line" | tr '\t' ' ') -> $(basename "$out")"
done

log "Wrote: $(basename "$SUMMARY")"
log "Now drawing figures"
Rscript "$SCRIPT_DIR/plot.R"
