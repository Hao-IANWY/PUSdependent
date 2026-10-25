#!/usr/bin/env bash
# data_prep.sh — trua3 Root：23 个跨组织 dependent 位点（cy-tRNA Ψ38/39/40）的逐重复测序深度
#
# 位点集合 = ../trua3_all_tissues.full.dependent.tsv 的 Root 行中 fix_pos ∈ {38,39,40}
# 的 cy 位点，即 *.full.* 分布图 Roots 组展示的那批位点。
# 深度来源（与 call_dependent_methylKit.R 进入检验的 WT/突变株 Treat 重复一一对应）：
#   WT      : WT-Root_sites.Treat_avg_ratio.tsv 的 Treat-WT-Root.rep{1,2,3}_depth
#   trua3_N : trua3_N-Root.methylkit.samples/trua3_N-Root.rep{1,2,3}.methyl.tsv 的 coverage
# 另输出每个样本在全部 439 个候选位点上的中位深度（证明覆盖度失衡是文库级的，
# 不限于这 23 个位点）。
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TRUA3_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../../../.." && pwd)"
SRC_DIR="$REPO_ROOT/tRNA/STAT/dependent_methylKit/ncRNA/Root/trua3"

FULL_TSV="$TRUA3_DIR/trua3_all_tissues.full.dependent.tsv"
OUT_TSV="$SCRIPT_DIR/trua3_root_union23_replicate_depth.tsv"
MEDIAN_TSV="$SCRIPT_DIR/trua3_root_candidate_median_depth.tsv"

readonly -A LINE_NAME=( ["1"]="trua3_1" ["2"]="trua3_2" )

[[ -s "$FULL_TSV" ]] || { echo "Missing full table: $FULL_TSV（先跑 ../data_prep.sh）" >&2; exit 1; }
[[ -s "$SRC_DIR/WT-Root_sites.Treat_avg_ratio.tsv" ]] || { echo "Missing WT candidate table: $SRC_DIR" >&2; exit 1; }

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

# ---------------------------------------------------------------------------
# 位点表：chr pos strand fix_pos（23 行）
# ---------------------------------------------------------------------------
SITE_TSV="$TMP_DIR/sites.tsv"
awk -F'\t' 'NR>1 && $11=="Root" && $3>=38 && $3<=40 {print $1"\t"$2"\t"$4"\t"$3}' \
  "$FULL_TSV" | sort -k1,1 -k2,2n > "$SITE_TSV"
n_sites="$(wc -l < "$SITE_TSV")"
[[ "$n_sites" -eq 23 ]] || { echo "Expected 23 sites, got $n_sites" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 长表：site fix_pos genotype rep depth（23 × 3 基因型 × 3 重复 = 207 行）
# ---------------------------------------------------------------------------
{
  printf 'site\tfix_pos\tgenotype\trep\tdepth\n'
  # WT 3 个 Treat 重复
  awk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] { site[$1"|"$2"|"$3] = $4; next }
    FNR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/\r$/, "", h); col[h] = i }; next }
    (($1"|"$2"|"$3) in site) {
      for (r = 1; r <= 3; r++)
        print $1 "_" $2, site[$1"|"$2"|"$3], "WT", "rep" r, $(col["Treat-WT-Root.rep" r "_depth"]) + 0
    }' "$SITE_TSV" "$SRC_DIR/WT-Root_sites.Treat_avg_ratio.tsv"
  # 两个株系各 3 个 Treat 重复
  for n in 1 2; do
    g="${LINE_NAME[$n]}"
    for r in 1 2 3; do
      f="$SRC_DIR/${g}-Root.methylkit.samples/${g}-Root.rep${r}.methyl.tsv"
      [[ -s "$f" ]] || { echo "Missing sample file: $f" >&2; exit 1; }
      awk -F'\t' -v OFS='\t' -v g="$g" -v rep="rep$r" '
        FILENAME == ARGV[1] { site[$1"|"$2"|"$3] = $4; next }
        (($1"|"$2"|"$4) in site) { print $1 "_" $2, site[$1"|"$2"|"$4], g, rep, $5 + 0 }
      ' "$SITE_TSV" "$f"
    done
  done
} > "$OUT_TSV"

n_rows="$(($(wc -l < "$OUT_TSV") - 1))"
[[ "$n_rows" -eq 207 ]] || { echo "Expected 207 depth rows, got $n_rows" >&2; exit 1; }

# ---------------------------------------------------------------------------
# 各样本在全部候选位点上的中位深度（439 行输入 = 439 个候选）
# ---------------------------------------------------------------------------
{
  printf 'sample\trep\tmedian_depth\n'
  for r in 1 2 3; do
    awk -F'\t' -v c="Treat-WT-Root.rep${r}_depth" '
      FNR == 1 { for (i = 1; i <= NF; i++) { h = $i; sub(/\r$/, "", h); col[h] = i }; next }
      { print $(col[c]) + 0 }' "$SRC_DIR/WT-Root_sites.Treat_avg_ratio.tsv" \
      | sort -g | awk -v s="WT" -v rep="rep$r" \
        '{ a[NR] = $1 } END { printf "%s\t%s\t%.0f\n", s, rep, (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2 }'
  done
  for n in 1 2; do
    g="${LINE_NAME[$n]}"
    for r in 1 2 3; do
      cat "$SRC_DIR/${g}-Root.methylkit.samples/${g}-Root.rep${r}.methyl.tsv" \
        | cut -f5 | sort -g \
        | awk -v s="$g" -v rep="rep$r" \
          '{ a[NR] = $1 } END { printf "%s\t%s\t%.0f\n", s, rep, (NR % 2) ? a[(NR + 1) / 2] : (a[NR / 2] + a[NR / 2 + 1]) / 2 }'
    done
  done
} > "$MEDIAN_TSV"

echo "[data_prep] Generated: $OUT_TSV ($n_rows rows, $n_sites sites × 9 sample-columns)"
echo "[data_prep] Generated: $MEDIAN_TSV ($(($(wc -l < "$MEDIAN_TSV") - 1)) sample-replicates)"
