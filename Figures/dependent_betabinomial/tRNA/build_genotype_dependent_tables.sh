#!/usr/bin/env bash
# build_genotype_dependent_tables.sh — Assemble per-tissue and cross-tissue tables of
# writer-dependent Ψ sites on tRNAs for one trua mutant genotype.
#
# Inputs:
#   tRNA/STAT/dependent_betabinomial/ncRNA/<tissue>/<genotype>/ — beta-binomial dependent
#   tables; Sprinzl positions must already be attached by
#   BID-transcript.script/downstream/tRNA_Sprinzl_annoation/add_sprinzl_pos.py
#   (file name <g>-<tissue>.dependent.Sprinzl.tsv).
#   <g>_{1,2}-<tissue>.dependent-sites.tsv — per-allele Mut_Ratio, joined by
#   chr+pos+strand. Mut_Ratio must come from these files (mean of the unfiltered adjustGap
#   pileup over the 3 Treat replicates of that allele), not from the Fisher tables,
#   because Fisher filtering zeroes weak residual mutant signals. No annotation columns
#   are needed here.
#   WT-Root_sites.Treat_avg_ratio.tsv — Root candidate universe (verified to be the same
#   quantity as WT_Ratio in the per-tissue dependent tables).
#
# Outputs (<output_dir>/):
#   <genotype>_all_tissues.dependent.tsv — per-tissue long table.
#   <genotype>_union.dependent.tsv — cross-tissue union.
#   <genotype>_all_tissues.full.dependent.tsv — leaf/Flower rows kept as is; Root rows
#     rebuilt from the union of the other two tissues with Root-measured fractions;
#     union sites outside the Root candidate universe (not measurable in Root) dropped.
#
# Notes:
#   fix_pos is taken directly from the upstream Sprinzl_pos column, so no coordinate
#   conversion is done here and All-tRNA_rep-structure.csv is not needed. Region
#   annotation is computed upstream per copy, so copies whose length or region layout
#   differs from the family representative are still handled correctly.
#   The shared ncRNA reference also contains sn/snoRNAs; keep only chr matching ^tRNA-.
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
readonly DEPENDENT_DIR="$ROOT_DIR/tRNA/STAT/dependent_betabinomial/ncRNA"

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
# Per-tissue long concatenation: dependent table (Sprinzl_pos included) + per-allele Mut_Ratio + Tissue
# ---------------------------------------------------------------------------
TMP_LONG="$TMP_DIR/all-tissues.tsv"
printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tMut_Ratio\tDelta_Ratio(WT-Mut)\t%s\t%s\tTissue\n' \
  "${GENOTYPE}_1" "${GENOTYPE}_2" > "$TMP_LONG"

for tissue in "${TISSUES[@]}"; do
  final_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/$GENOTYPE-$tissue.dependent.Sprinzl.tsv"
  line1_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_1-$tissue.dependent-sites.tsv"
  line2_file="$DEPENDENT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_2-$tissue.dependent-sites.tsv"

  [[ -s "$final_file" ]] || error "Missing annotated dependent table: $final_file (run add_sprinzl_pos.py first)"
  [[ -s "$line1_file" ]] || error "Missing line-level table: $line1_file"
  [[ -s "$line2_file" ]] || error "Missing line-level table: $line2_file"

  awk -F'\t' -v OFS='\t' \
      -v tissue="$tissue" \
      -v trna_pat="$TRNA_CHR_PATTERN" \
      -v line1_file="$line1_file" \
      -v line2_file="$line2_file" '
  # Column name -> column index (the three inputs may differ in column order; always resolve from the header)
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
  # Third input: the dependent table (carries Sprinzl_pos)
  FNR == 1 { load_cols(); require_cols("chr,pos,Sprinzl_pos,strand,context,WT_Ratio,Mut_Ratio,Delta_Ratio(WT-Mut)"); next }
  {
      chr = $(col["chr"])
      if (chr !~ trna_pat) next                  # Drop sn/snoRNAs that share the ncRNA reference

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
# Cross-tissue union: keep the observation with the largest Delta per site; add Tissues / n_tissues / peak_tissue
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

# ---------------------------------------------------------------------------
# .full table: keep leaf/Flower rows unchanged; replace Root rows with the union of the
# dependent sites of the other two tissues, with all modification fractions measured in
# Root — WT from WT-Root_sites.Treat_avg_ratio (verified to be the same quantity as the
# WT_Ratio of the per-tissue dependent tables), alleles from their own Root
# dependent-sites Mut_Ratio; union sites outside the Root candidate universe (not
# measurable in Root) are dropped.
# ---------------------------------------------------------------------------
FULL_OUT="$OUT_DIR/${GENOTYPE}_all_tissues.full.dependent.tsv"

WT_ROOT_FILE="$DEPENDENT_DIR/Root/$GENOTYPE/WT-Root_sites.Treat_avg_ratio.tsv"
ROOT_LINE1_FILE="$DEPENDENT_DIR/Root/$GENOTYPE/${GENOTYPE}_1-Root.dependent-sites.tsv"
ROOT_LINE2_FILE="$DEPENDENT_DIR/Root/$GENOTYPE/${GENOTYPE}_2-Root.dependent-sites.tsv"
[[ -s "$WT_ROOT_FILE" ]] || error "Missing WT Root candidate table: $WT_ROOT_FILE"
[[ -s "$ROOT_LINE1_FILE" ]] || error "Missing line-level table: $ROOT_LINE1_FILE"
[[ -s "$ROOT_LINE2_FILE" ]] || error "Missing line-level table: $ROOT_LINE2_FILE"

{
  printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tMut_Ratio\tDelta_Ratio(WT-Mut)\t%s\t%s\tTissue\n' \
    "${GENOTYPE}_1" "${GENOTYPE}_2"
  awk -F'\t' -v OFS='\t' \
      -v genotype="$GENOTYPE" \
      -v wt_root_file="$WT_ROOT_FILE" \
      -v line1_file="$ROOT_LINE1_FILE" \
      -v line2_file="$ROOT_LINE2_FILE" '
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
  FILENAME == wt_root_file {
      if (FNR == 1) { load_cols(); require_cols("chr,pos,strand,Treat_avg_ratio"); next }
      wt_root[$1 SUBSEP $2 SUBSEP $3] = $(col["Treat_avg_ratio"])
      next
  }
  FILENAME == line1_file {
      if (FNR == 1) { load_cols(); require_cols("chr,pos,strand,Mut_Ratio"); next }
      mut1_root[$1 SUBSEP $2 SUBSEP $3] = $(col["Mut_Ratio"])
      next
  }
  FILENAME == line2_file {
      if (FNR == 1) { load_cols(); require_cols("chr,pos,strand,Mut_Ratio"); next }
      mut2_root[$1 SUBSEP $2 SUBSEP $3] = $(col["Mut_Ratio"])
      next
  }
  # Last input: the all_tissues long table
  FNR == 1 { next }
  $NF == "Root" { next }                     # Discard the original Root rows; they are rebuilt from the union
  {
      print                                  # leaf/Flower rows pass through unchanged
      skey = $1 SUBSEP $2 SUBSEP $4
      if (!(skey in union_fix)) {            # Deduplicate across tissues; record first fix_pos/context seen
          union_order[++n_union] = skey
          union_fix[skey] = $3
          union_ctx[skey] = $5
      }
      next
  }
  END {
      n_dropped = 0
      for (i = 1; i <= n_union; i++) {
          skey = union_order[i]
          if (!(skey in wt_root)) { n_dropped++; continue }
          split(skey, f, SUBSEP)
          mut1 = (skey in mut1_root) ? mut1_root[skey] : 0
          mut2 = (skey in mut2_root) ? mut2_root[skey] : 0
          pooled_mut = (mut1 + mut2) / 2
          print f[1], f[2], union_fix[skey], f[3], union_ctx[skey], \
                wt_root[skey], pooled_mut, wt_root[skey] - pooled_mut, \
                mut1, mut2, "Root"
      }
      printf("[data_prep] %s: .full Root rows = %d union sites (%d dropped: no Root measurement)\n", \
             genotype, n_union - n_dropped, n_dropped) > "/dev/stderr"
  }' "$WT_ROOT_FILE" "$ROOT_LINE1_FILE" "$ROOT_LINE2_FILE" "$TMP_LONG"
} > "$FULL_OUT"

log "$GENOTYPE: $(( $(wc -l < "$LONG_OUT") - 1 )) tissue-resolved rows, $(( $(wc -l < "$UNION_OUT") - 1 )) union sites, $(( $(wc -l < "$FULL_OUT") - 1 )) .full rows"
