#!/usr/bin/env bash
# data.sh — 用 call_dependent_betabinomial.R 做 Root/Leaf(参照) vs Flower(mutant) 的
# 组织间 tRNA 修饰比例比较，跑两次取并集得 Flower-low 位点，再交 merge_flower_low.R /
# plot.R 出表出图。口径与细节见 README.md。
#
# 用法:
#   bash data.sh              # 全流程（两次 bb 调用并行 + 合并 + 绘图）
#   bash data.sh --stats-only # 只跑统计与合并
#   bash data.sh --plot-only  # 只重画图（输入表已在）
set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# SCRIPT_DIR = <repo>/Figures/wildtype_landscape/tissue_ratio/tRNA → 4 层回到仓库根
readonly REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
readonly RSCRIPT_BIN="/home/HL/miniforge3/envs/analysis/bin/Rscript"
readonly CALL_SCRIPT="${REPO_ROOT}/BID-transcript.script/downstream/call_dependent_betabinomial.R"
readonly FISHER_DIR="${REPO_ROOT}/tRNA/STAT/fisher/ncRNA"
readonly ADJUST_DIR="${REPO_ROOT}/tRNA/STAT/adjustGap/ncRNA"
readonly BB_DIR="${SCRIPT_DIR}/betabinomial"

# 组织比较口径（见 README）：q < 0.05 且 delta >= 0.05 才算 Flower-low
# （老板 2026-11-03 拍板：与其他 dependent 位点分析的效应量门槛一致）。
readonly MIN_DELTA="0.05"
readonly ALPHA="0.05"
readonly MIN_DEPTH="10"
readonly CORES="4"

log() { printf '[data] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
error() { log "ERROR: $*"; exit 1; }

[[ -x "$RSCRIPT_BIN" ]] || error "Rscript not found: ${RSCRIPT_BIN}"
[[ -f "$CALL_SCRIPT" ]] || error "bb script not found: ${CALL_SCRIPT}"

# WT 角色 --tsv：Fisher 表 Pass==1 且 ^tRNA- 的行（BH 家族 = tRNA 候选池，
# sn/snoRNA 不参与校正，口径与 dependent bb 一致）。
prepare_wt_tsv() {
    local fisher_tsv="$1" output_tsv="$2"
    awk 'BEGIN{FS=OFS="\t"} NR == 1 || ($NF == 1 && $1 ~ /^tRNA-/)' "$fisher_tsv" > "$output_tsv"
}

# 跑一次组织比较：ref_tissue 为参照组织（WT 角色），Flower 恒为 mutant 角色。
# genotype2 与 genotype1 传同一组 Flower 文件（脚本双株系接口的占位，AND 退化为单比较）。
run_comparison() {
    local ref_tissue="$1"          # Root | Leaf
    local outdir="${BB_DIR}/${ref_tissue,,}_vs_flower"
    local genotype="Flower_vs_${ref_tissue}"
    local wt_tsv="${outdir}/WT-${ref_tissue}_sites.Treat_avg_ratio.tsv"

    mkdir -p "$outdir"
    prepare_wt_tsv "${FISHER_DIR}/WT_${ref_tissue}_sites.Fisher.tsv" "$wt_tsv"
    log "(${ref_tissue,,}_vs_flower) WT candidates: $(($(wc -l < "$wt_tsv") - 1)) Pass tRNA sites"

    local -a treat_files=("${ADJUST_DIR}/Treat-WT-Flower.rep1.adjusted.tsv" \
                          "${ADJUST_DIR}/Treat-WT-Flower.rep2.adjusted.tsv" \
                          "${ADJUST_DIR}/Treat-WT-Flower.rep3.adjusted.tsv")
    local -a mock_files=("${ADJUST_DIR}/Mock-WT-Flower.rep1.adjusted.tsv" \
                          "${ADJUST_DIR}/Mock-WT-Flower.rep2.adjusted.tsv" \
                          "${ADJUST_DIR}/Mock-WT-Flower.rep3.adjusted.tsv")
    local file
    for file in "${treat_files[@]}" "${mock_files[@]}"; do
        [[ -s "$file" ]] || error "File not found: ${file}"
    done
    local treat_csv mock_csv
    treat_csv="$(IFS=,; printf '%s' "${treat_files[*]}")"
    mock_csv="$(IFS=,; printf '%s' "${mock_files[*]}")"

    "$RSCRIPT_BIN" "$CALL_SCRIPT" \
        --tsv "$wt_tsv" \
        --key-cols chr,pos,strand \
        --ratio-col Treat_avg_ratio \
        --significance q \
        --site-test lrt \
        --alpha "$ALPHA" \
        --min-depth "$MIN_DEPTH" \
        --min-delta "$MIN_DELTA" \
        --cores "$CORES" \
        --genotype1 "$genotype" \
        --genotype1-treat-files "$treat_csv" \
        --genotype1-mock-files "$mock_csv" \
        --genotype2 "$genotype" \
        --genotype2-treat-files "$treat_csv" \
        --genotype2-mock-files "$mock_csv" \
        --outdir "$outdir" \
        --output-prefix "${ref_tissue,,}_vs_flower" \
        > "${outdir}/run.log" 2>&1
}

MODE="all"
case "${1:-all}" in
    --stats-only) MODE="stats" ;;
    --plot-only)  MODE="plot" ;;
    all)          MODE="all" ;;
    *)            error "Unknown option: $1 (use --stats-only / --plot-only)" ;;
esac

if [[ "$MODE" != "plot" ]]; then
    run_comparison Root &
    pid_root=$!
    run_comparison Leaf &
    pid_leaf=$!
    fail=0
    wait "$pid_root" || { log "FAILED root_vs_flower; see ${BB_DIR}/root_vs_flower/run.log"; fail=1; }
    wait "$pid_leaf" || { log "FAILED leaf_vs_flower; see ${BB_DIR}/leaf_vs_flower/run.log"; fail=1; }
    [[ "$fail" -eq 0 ]] || error "beta-binomial comparison task(s) failed"
    log "bb comparisons done; merging Flower-low sites"
    "$RSCRIPT_BIN" "${SCRIPT_DIR}/merge_flower_low.R"
fi

if [[ "$MODE" != "stats" ]]; then
    log "drawing heatmap"
    "$RSCRIPT_BIN" "${SCRIPT_DIR}/plot.R"
fi

log "Done. Outputs are under ${SCRIPT_DIR}"
