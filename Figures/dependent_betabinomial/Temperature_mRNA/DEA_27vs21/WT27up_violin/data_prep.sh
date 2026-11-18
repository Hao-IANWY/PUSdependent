#!/usr/bin/env bash
# ============================================================================
# data_prep.sh — WT 27°C 上调基因在 {WT, trua1b_1, trua1b_2} 的 27°C 表达量
#
# 基因集：WT_27C_vs_21C_DESeq2.csv 中 padj < 0.05 且 log2FoldChange >= 1
#         的 Upregulated 基因（891 个；与上级火山图 LOGFC_FILTER=1 /
#         PADJ_FILTER=0.05 的默认阈值一致，可用环境变量覆盖）
# 表达量：各基因型 27°C 的 DESeq2 归一化 count（rep3/4/5）
#
# 缺失基因还原：各突变体在 DESeq2 前做了"组内 6 样品 count 全 > 0"过滤，
#   少数 WT 上调基因（trua1b_1 缺 37、trua1b_2 缺 41 个，均为混合零）不在
#   其结果 CSV 中。DESeq2 归一化定义 norm = raw / sizeFactor，故对任一
#   CSV 内基因有 sizeFactor = raw / norm；据此按样品反推 size factor
#   （全 CSV 基因上取均值，规避浮点噪声），再把缺失基因按 norm = raw / sf
#   从 count_27vs21.tsv 精确还原，保证 891 个基因在三组中完整参与展示。
#
# 输出：WT27up_expr_27C.tsv（长表：genotype / sample / feature_id / norm_count）
#       891 基因 × 3 样品 × 3 基因型 = 8019 行
#
# 用法：
#   bash data_prep.sh
#   LOGFC_FILTER=1 PADJ_FILTER=0.05 bash data_prep.sh   # 覆盖阈值
# ============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="${SCRIPT_DIR}/.."                       # 上级 DEA_27vs21 目录
COUNT_TSV="${BASE_DIR}/count_27vs21.tsv"
LOGFC_FILTER="${LOGFC_FILTER:-1}"
PADJ_FILTER="${PADJ_FILTER:-0.05}"
OUT_TSV="${SCRIPT_DIR}/WT27up_expr_27C.tsv"
GENOTYPES=(WT trua1b_1 trua1b_2)

[ -f "$COUNT_TSV" ] || { echo "[data_prep] ERROR: count matrix not found: $COUNT_TSV" >&2; exit 1; }

# ---- 1. WT 上调基因列表（padj < PADJ_FILTER 且 log2FC >= LOGFC_FILTER）----
GENE_LIST="$(mktemp)"
trap 'rm -f "$GENE_LIST"' EXIT
awk -F',' -v lfc="$LOGFC_FILTER" -v padj="$PADJ_FILTER" '
  NR > 1 && $7 != "NA" && $7 != "" && $7+0 < padj && $3+0 >= lfc { print $1 }
' "${BASE_DIR}/WT/WT_27C_vs_21C_DESeq2.csv" | sort -u > "$GENE_LIST"

N_GENES=$(wc -l < "$GENE_LIST")
[ "$N_GENES" -gt 0 ] || { echo "[data_prep] ERROR: no up genes selected" >&2; exit 1; }

# ---- 2. 逐基因型提取 27°C 归一化 count（按表头列名定位，缺列即报错）----
printf "genotype\tsample\tfeature_id\tnorm_count\n" > "$OUT_TSV"
for g in "${GENOTYPES[@]}"; do
  csv="${BASE_DIR}/${g}/${g}_27C_vs_21C_DESeq2.csv"
  [ -f "$csv" ] || { echo "[data_prep] ERROR: missing $csv" >&2; exit 1; }
  awk -F',' -v OFS='\t' \
      -v counts="$COUNT_TSV" -v genes="$GENE_LIST" -v geno="$g" '
    # --- 预载：上调基因集合 + count 矩阵中本基因型 27°C 三列的原始 count ---
    BEGIN {
      while ((getline line < genes) > 0) keep[line] = 1
      close(genes)
      ncol = 0
      while ((getline line < counts) > 0) {
        m = split(line, a, "\t")
        if (a[1] == "feature_id") {
          for (i = 1; i <= m; i++) cidx[a[i]] = i
          for (r = 3; r <= 5; r++) {
            col = "Mock-27-" geno ".rep" r
            if (!(col in cidx)) {
              printf "[data_prep] ERROR: count matrix missing column %s\n", col > "/dev/stderr"
              failed = 1; exit 1
            }
            cpos[ncol] = cidx[col]; cname[ncol] = col; ncol++
          }
          continue
        }
        for (j = 0; j < ncol; j++) raw[a[1], j] = a[cpos[j]] + 0
      }
      close(counts)
    }
    # --- CSV 表头：定位本基因型 27°C 三列 ---
    NR == 1 {
      for (i = 1; i <= NF; i++) vidx[$i] = i
      for (j = 0; j < ncol; j++) {
        if (!(cname[j] in vidx)) {
          printf "[data_prep] ERROR: %s missing column %s\n", geno, cname[j] > "/dev/stderr"
          failed = 1; exit 1
        }
        vpos[j] = vidx[cname[j]]
      }
      next
    }
    # --- 数据行：CSV 内基因累计 size factor 比值；上调基因记录归一化值 ---
    {
      id = $1
      for (j = 0; j < ncol; j++) {
        v = $vpos[j]
        if (v != "NA" && v != "" && (v + 0) == v) {
          rv = raw[id, j]
          if (rv > 0) { sfsum[j] += rv / (v + 0); sfn[j]++ }
          if (id in keep) { upnorm[id, j] = v + 0; have[id, j] = 1 }
        }
      }
    }
    # --- 输出：CSV 内基因用原值，缺失基因按 norm = raw / sf 还原 ---
    END {
      if (failed) exit 1
      for (j = 0; j < ncol; j++) {
        if (sfn[j] == 0) {
          printf "[data_prep] ERROR: cannot infer size factor for %s\n", cname[j] > "/dev/stderr"
          exit 1
        }
        sfv[j] = sfsum[j] / sfn[j]
      }
      while ((getline id < genes) > 0) {
        for (j = 0; j < ncol; j++) {
          if (have[id, j]) v = upnorm[id, j]
          else             v = raw[id, j] / sfv[j]
          printf "%s\t%s\t%s\t%.6f\n", geno, cname[j], id, v
        }
      }
      close(genes)
    }
  ' "$csv" >> "$OUT_TSV"
done

N_ROWS=$(( $(wc -l < "$OUT_TSV") - 1 ))
echo "[data_prep] WT 27C-up genes: $N_GENES; exported rows: $N_ROWS -> $OUT_TSV"
