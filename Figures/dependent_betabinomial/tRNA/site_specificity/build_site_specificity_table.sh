#!/usr/bin/env bash
# build_site_specificity_table.sh — Build the data table for the site-specificity analysis.
#
# Question: how many representative tRNAs carry modification at a given canonical
# (Sprinzl/fixed) position in WT, and what modification fractions do those tRNAs show in
# a chosen set of target mutants (two alleles each)?
#
# Inclusion criteria (all four must hold):
#   1. The site lies in the WT Fisher-positive universe (WT-<tissue>_sites.Treat_avg_ratio.tsv;
#      every row of that table has Pass=1, i.e. the CMC effect was detected in all three
#      WT Treat replicates = modified in WT);
#   2. Sprinzl_pos == the requested position (native pos differ; 26/27/28/29 may all map
#      to 27);
#   3. chr is the family representative tRNA (the highest mean_CPM member of each
#      family);
#   4. The compartment prefix is on the compartments allow-list (empty = all
#      compartments).
#
# Fraction criteria: WT uses Treat_avg_ratio; each allele uses the Mut_Ratio of
#   <g>_<line>-<tissue>.dependent-sites.tsv (the mean over that allele's Treat
#   replicates, unfiltered), matching the per-genotype figure set; a site missing from
#   an allele table stays NA (never treated as 0% modification), and the per-line
#   is_dependent is carried alongside.
#   The WT table is byte-identical across genotype directories, so the first genotype
#   directory's copy is used.
#
# The shared ncRNA reference mixes in sn/snoRNAs; drop them with ^tRNA- (same criteria
# as build_genotype_dependent_tables.sh).
#
# Usage: bash build_site_specificity_table.sh <tissue> <genotypes> <sprinzl_pos> <output_dir> [compartments]
#   tissue is one of Root|leaf|Flower; genotypes is a comma list (e.g. trua1a,trua5);
#   sprinzl_pos e.g. 27; compartments is a comma list of chr-prefix compartments
#   (Nu/Mt/Pt); empty = all.
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
# This directory sits one level deeper than the genotype directory (the layer holding
# build_genotype_dependent_tables.sh), so step up one extra level
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

# Two allele tables per genotype; also verify the WT table in each genotype directory
# is byte-identical to the copy in use
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

# The raw WT-positive universe lacks Sprinzl_pos (annotation was only run on the 21
# final dependent tables); an annotated intermediate copy is produced here in the shared
# layer (this script's directory): it is reused by all position folders of the same
# tissue, and the tRNA/STAT/ raw data stays untouched.
readonly WT_SPRINZL_FILE="$SCRIPT_DIR/WT-${TISSUE}_sites.Sprinzl.tsv"
readonly STEM_TISSUE="$(printf '%s' "$TISSUE" | tr '[:upper:]' '[:lower:]')"
readonly OUT_TSV="$OUT_DIR/pos${SPRINZL_POS}_${STEM_TISSUE}_site_table.tsv"

python3 "$SPRINZL_SCRIPT" -i "$WT_FILE" -o "$WT_SPRINZL_FILE" >&2

# Header: two allele-value columns + two per-line dependent columns, expanded in genotype order
{
  printf 'chr\tpos\tfix_pos\tstrand\tcontext\tWT_Ratio\tTreat_avg_depth'
  for g in "${GENOTYPE_LIST[@]}"; do
    printf '\t%s_1\t%s_2\t%s_1_dependent\t%s_2_dependent' "$g" "$g" "$g" "$g"
  done
  printf '\n'
} > "$OUT_TSV"

# The allele-table path list is passed to awk; the rank of FILENAME in the list maps to
# (genotype, allele)
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
# Column name -> column index (input column order is not guaranteed; always resolve from the header)
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
    # Compartment restriction (empty = no filtering): chr prefixes tRNA-Nu-/tRNA-Mt-/tRNA-Pt- map to Nu/Mt/Pt
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
# Allele tables: rank of FILENAME in line_file_arg -> (genotype index gi, allele ln)
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
            # Column order matches the header: the two allele values first, then the two per-line dependent flags
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
