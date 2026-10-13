#!/usr/bin/env bash
# data.sh — 把 WT 各组织的 tRNA positive site 投影到 canonical tRNA 结构布局,
#           供 plot.R 画代表性核 tRNA 热图。
#
# 输入:
#   tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}_sites.positive.tsv
#   Nu-tRNA_rep-structure.csv          46 个核代表的固定行序(来自 archive,见 README)
#   Nu-tRNA_rep_structure_regions.tsv  每个代表的区段序列与坐标(含 CCA 尾的 1-based 参考坐标)
#   /home/HL/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa  仅用于识别长度与代表不符的拷贝
#
# 输出:
#   representative_nuclear_trna_mapped_sites_{root,leaf,flower}.tsv
#       tissue, Rep_ID, chr, pos, strand, ratio, depth,
#       region, region_pos, fixed_pos, matrix_col_id, canonical_pos
#   trna_landscape_summary.tsv
#
# 投影规则(照 archive 论文版实现):
#   位点 pos 落在某区段 [start_pos, end_pos] 内 -> region_pos = pos - start_pos + 1
#   再按区段固定长度重标定到 canonical 位:
#     round((region_pos - 1) * (fixed_length - 1) / (region_length - 1)) + 1
#   region_length 或 fixed_length 为 1 时取 1。canonical_pos 为 0-76 的展示位号。
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
readonly INPUT_DIR="$REPO_ROOT/tRNA/STAT/fisher/ncRNA"
readonly REGIONS="$SCRIPT_DIR/Nu-tRNA_rep_structure_regions.tsv"
readonly COLUMNS=$'tissue\tRep_ID\tchr\tpos\tstrand\tratio\tdepth\tregion\tregion_pos\tfixed_pos\tmatrix_col_id\tcanonical_pos'

log() { printf '[data_prep] %s\n' "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

# 参考 fasta:与 tRNA/config.yaml 的 ncRNA 参考一致
readonly -a FASTA_CANDIDATES=(
    /home/HL/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa
    /home/haoli/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa
)
FASTA=""
for f in "${FASTA_CANDIDATES[@]}"; do
    [[ -f "$f" ]] && { FASTA="$f"; break; }
done
[[ -n "$FASTA" ]] || error "tRNA reference fasta not found in: ${FASTA_CANDIDATES[*]}"
[[ -s "$REGIONS" ]] || error "Missing structure table: $REGIONS"

readonly SUMMARY="$SCRIPT_DIR/trna_landscape_summary.tsv"
readonly COUNTS="$SCRIPT_DIR/.counts.tmp"
trap 'rm -f "$COUNTS"' EXIT

printf 'tissue\tsites_in\texcluded_len_mismatch\tunmapped\tsites_mapped\treps_covered\tcells\n' > "$SUMMARY"

for tissue in Root Leaf Flower; do
    input="$INPUT_DIR/WT_${tissue}_sites.positive.tsv"
    [[ -s "$input" ]] || error "Missing input: $input"
    prefix="$(printf '%s' "$tissue" | tr '[:upper:]' '[:lower:]')"
    out="$SCRIPT_DIR/representative_nuclear_trna_mapped_sites_${prefix}.tsv"

    printf '%s\n' "$COLUMNS" > "$out"
    awk -F'\t' -v OFS='\t' -v tissue="$tissue" -v regions="$REGIONS" -v fasta="$FASTA" -v counts="$COUNTS" '
    # canonical 布局:区段顺序、固定长度、在 0-76 展示位中的起始偏移
    BEGIN {
        split("p1 Accstem_5p p8-9 Dstem_5p Dloop Dstem_3p p26 Acstem_5p Acloop Acstem_3p " \
              "Vregion Tstem_5p Tloop Tstem_3p Accstem_3p p73 CCA", order, " ")
        split("1 7 2 4 8 4 1 5 7 5 5 5 7 5 7 1 3", lens, " ")
        off = 0
        for (i = 1; i <= 17; i++) { FIXLEN[order[i]] = lens[i]; OFFSET[order[i]] = off; off += lens[i] }
        ncol = off
    }
    # R 的 round() 是银行家舍入(0.5 取偶),这里必须一致,否则 .5 边界会串列
    function round_like_r(x,   fl) {
        fl = int(x)
        if (x - fl == 0.5) return (fl % 2 == 0) ? fl : fl + 1
        return int(x + 0.5)
    }
    FILENAME == regions {
        if (FNR > 1 && $5 + 0 > 0) {
            n++; rep = $1
            R_rep[n] = rep; R_region[n] = $3; R_len[n] = $5 + 0
            R_st[n] = $6 + 0; R_en[n] = $7 + 0
            if (R_en[n] > span[rep]) span[rep] = R_en[n]
        }
        next
    }
    FILENAME == fasta {
        if (/^>/) {
            name = substr($0, 2); sub(/[ \t].*$/, "", name); seq = ""
        } else {
            seq = seq $0
            if (name ~ /^tRNA-Nu-/) copy_len[name] = length(seq)
        }
        next
    }
    # 第三个输入:positive sites(有表头)
    FNR == 1 { next }
    {
        sites_in++
        chr = $1
        if (chr !~ /^tRNA-Nu-/) next          # 只画核 tRNA
        rep = chr; sub(/^tRNA-Nu-/, "Nu-", rep); sub(/-[0-9]+$/, "", rep)
        if (rep in span && copy_len[chr] != "" && copy_len[chr] != span[rep]) {
            excluded++; next                  # 拷贝序列长度与代表不符,坐标整体位移
        }
        pos = $2 + 0; ratio = $4 + 0; depth = $5 + 0
        hit = 0; sel = 0
        for (i = 1; i <= n; i++) {
            if (R_rep[i] == rep && pos >= R_st[i] && pos <= R_en[i]) { hit++; sel = i }
        }
        if (hit != 1) { unmapped++; next }
        region = R_region[sel]
        region_pos = pos - R_st[sel] + 1
        fixed_len = FIXLEN[region]
        if (fixed_len == 1 || R_len[sel] == 1) fixed_pos = 1
        else fixed_pos = round_like_r((region_pos - 1) * (fixed_len - 1) / (R_len[sel] - 1)) + 1
        colid = region "_" fixed_pos

        print tissue, rep, chr, pos, $3, sprintf("%.10g", ratio), sprintf("%.10g", depth),
              region, region_pos, fixed_pos, colid, OFFSET[region] + fixed_pos - 1
        mapped++
        seen_rep[rep] = 1
        seen_cell[rep SUBSEP colid] = 1
    }
    END {
        printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", tissue, sites_in + 0, excluded + 0,
               unmapped + 0, mapped + 0, length(seen_rep), length(seen_cell) > counts
    }
    ' "$REGIONS" "$FASTA" "$input" | sort -t$'\t' -k2,2 -k12,12n -k4,4n >> "$out"

    counts_line="$(cat "$COUNTS")"
    printf '%s\n' "$counts_line" >> "$SUMMARY"
    log "$(printf '%s' "$counts_line" | tr '\t' ' ') -> $(basename "$out")"
done

log "Wrote: $(basename "$SUMMARY")"
log "Now drawing figures"
Rscript "$SCRIPT_DIR/plot.R"
