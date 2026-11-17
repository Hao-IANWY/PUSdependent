#!/usr/bin/env bash
# ============================================================================
# data_prep.sh — 27°C vs 21°C 差异分析数据整备（WT / trua1b_1 / trua1b_2，Mock）
#
# 取值：Mock 侧 rep3/rep4/rep5（rep1/rep2 与 rep3 属混批设计，改用 rep3-5
# 同批组合；上游为含 rep4/5 的 Mock.primary_counts.tsv）
#
# 职责（只整备数据，不画图）：
#   1. 从 Mock.primary_counts.tsv 提取三个基因型 × 21/27°C × rep3-5 共 18 列
#      → count_27vs21.tsv（行=基因，第一列=feature_id）
#   2. 写样品元数据 samples.tsv（sample / genotype / temperature）
#
# 用法：
#   bash data_prep.sh [count源文件]
#   默认源文件经仓库根的 Temperature 符号链接指向
#   <repo>/Temperature/count/mRNA/Mock.primary_counts.tsv
#   （本目录位于 Figures/dependent_betabinomial/Temperature_mRNA/ 下，
#    相对路径 ../../../../Temperature/... = 仓库根；整体挪动目录树时只需改这一处）
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC="$(readlink -f "${1:-${SCRIPT_DIR}/../../../../Temperature/count/mRNA/Mock.primary_counts.tsv}")"
OUT_COUNT="${SCRIPT_DIR}/count_27vs21.tsv"
OUT_META="${SCRIPT_DIR}/samples.tsv"

GENOTYPES=(WT trua1b_1 trua1b_2)   # 面板顺序
TEMPERATURES=(21 27)               # 21=对照, 27=处理
REPS=(rep3 rep4 rep5)              # 同批重复组合

[ -f "$SRC" ] || { echo "[data_prep] ERROR: source not found: $SRC" >&2; exit 1; }

# ---- 组装目标列（基因型为主序，温度为次序）----
COLS="feature_id"
for g in "${GENOTYPES[@]}"; do
  for t in "${TEMPERATURES[@]}"; do
    for r in "${REPS[@]}"; do
      COLS="${COLS} Mock-${t}-${g}.${r}"
    done
  done
done

# ---- 按表头列名提取；缺列即报错 ----
awk -v COLS="$COLS" -F'\t' '
  BEGIN { n = split(COLS, want, " ") }
  NR == 1 {
    for (i = 1; i <= NF; i++) idx[$i] = i
    miss = ""
    for (j = 1; j <= n; j++) if (!(want[j] in idx)) miss = miss "  " want[j] "\n"
    if (miss != "") {
      printf "[data_prep] ERROR: columns missing from source:\n%s", miss > "/dev/stderr"
      exit 1
    }
    for (j = 1; j <= n; j++) printf "%s%s", (j > 1 ? "\t" : ""), want[j]
    printf "\n"
    next
  }
  {
    for (j = 1; j <= n; j++) printf "%s%s", (j > 1 ? "\t" : ""), $idx[want[j]]
    printf "\n"
  }' "$SRC" > "$OUT_COUNT"

# ---- 样品元数据 ----
{
  printf "sample\tgenotype\ttemperature\n"
  for g in "${GENOTYPES[@]}"; do
    for t in "${TEMPERATURES[@]}"; do
      for r in "${REPS[@]}"; do
        printf "Mock-%s-%s.%s\t%s\t%s\n" "$t" "$g" "$r" "$g" "$t"
      done
    done
  done
} > "$OUT_META"

n_genes=$(( $(wc -l < "$OUT_COUNT") - 1 ))
n_samples=$(( $(head -1 "$OUT_COUNT" | awk -F'\t' '{print NF}') - 1 ))
echo "[data_prep] Count matrix: $OUT_COUNT  ($n_genes genes x $n_samples samples)"
echo "[data_prep] Sample meta : $OUT_META  ($(( $(wc -l < "$OUT_META") - 1 )) samples)"
