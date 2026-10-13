#!/usr/bin/env bash
# data_prep.sh — Build the cross-tissue site union and six-group ratio table for rRNA trua1a
#
# Inputs:
#   polyARNA/STAT/dependent_methylKit/rRNA/<tissue>/trua1a/trua1a-<tissue>.dependent.tsv
#   polyARNA/STAT/dependent_methylKit/rRNA/<tissue>/trua1a/WT-<tissue>_sites.Treat_avg_ratio.tsv
#   polyARNA/STAT/dependent_methylKit/rRNA/<tissue>/trua1a/trua1a_{1,2}-<tissue>.dependent-sites.tsv
#
# Outputs:
#   tissue_specificity_union_sites.tsv
#   tissue_specificity_ratio_summary.tsv
set -euo pipefail
IFS=$'\n\t'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
INPUT_DIR="$ROOT_DIR/polyARNA/STAT/dependent_methylKit/rRNA"
GENOTYPE="trua1a"
UNION_TSV="$SCRIPT_DIR/tissue_specificity_union_sites.tsv"
RATIO_TSV="$SCRIPT_DIR/tissue_specificity_ratio_summary.tsv"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

tissues=("Root" "leaf" "Flower")
tissue_labels=("Root" "Leaf" "Flower")

require_file() {
    local file="$1"
    [[ -s "$file" ]] || {
        echo "[data_prep] Missing or empty input file: $file" >&2
        exit 1
    }
}

# ---------------------------------------------------------------------------
# 1. Site union: all trua1a-dependent rRNA sites across Root, leaf, and Flower.
# ---------------------------------------------------------------------------
dependent_files=()
for tissue in "${tissues[@]}"; do
    dependent_file="$INPUT_DIR/$tissue/$GENOTYPE/${GENOTYPE}-${tissue}.dependent.tsv"
    require_file "$dependent_file"
    dependent_files+=("$dependent_file")
done

# A single awk process keeps one in-memory site dictionary, so sites shared by
# more than one tissue are collapsed into one row and their source tissues are
# retained in deterministic Root,Leaf,Flower order.
awk -F '\t' -v OFS='\t' -v genotype="$GENOTYPE" -v tissue_arg="Root,Leaf,Flower" '
    BEGIN {
        split(tissue_arg, tissue_labels, ",")
        split("chr,pos,strand,context", required_names, ",")
        for (i in required_names) required[required_names[i]] = 1
    }
    FNR == 1 {
        file_index++
        delete col
        for (i = 1; i <= NF; i++) {
            sub(/\r$/, "", $i)
            col[$i] = i
        }
        for (name in required) {
            if (!(name in col)) {
                printf("Required column missing in %s: %s\n", FILENAME, name) > "/dev/stderr"
                exit 1
            }
        }
        next
    }
    {
        tissue_label = tissue_labels[file_index]
        chr = $(col["chr"]); pos = $(col["pos"]); strand = $(col["strand"])
        context = $(col["context"])
        sub(/\r$/, "", context)
        key = chr "|" pos "|" strand
        if (!(key in seen)) {
            seen[key] = 1
            chr_map[key] = chr; pos_map[key] = pos; strand_map[key] = strand
            context_map[key] = context; tissues_map[key] = ""; n_map[key] = 0
        } else if (context_map[key] != context) {
            printf("Context mismatch for %s: %s vs %s\n", key, context_map[key], context) > "/dev/stderr"
            exit 1
        }
        if (!((key, tissue_label) in tissue_seen)) {
            tissue_seen[key, tissue_label] = 1
            n_map[key]++
            tissues_map[key] = (tissues_map[key] == "" ? tissue_label : tissues_map[key] "," tissue_label)
        }
    }
    END {
        for (key in seen) {
            print genotype, chr_map[key], pos_map[key], strand_map[key], context_map[key], key, tissues_map[key], n_map[key]
        }
    }
' "${dependent_files[@]}" > "$TMP_DIR/union.body.tsv"

{
    printf 'genotype\tchr\tpos\tstrand\tcontext\tsite_id\tsource_tissues\tn_tissues\n'
    sort -t "$(printf '\t')" -k2,2 -k3,3n -k4,4 "$TMP_DIR/union.body.tsv"
} > "$UNION_TSV"


# ---------------------------------------------------------------------------
# 2. Six-group ratio table. The same site union is used in every group.
#    WT: Pass == 1 and Treat_avg_ratio (three WT replicates, already averaged).
#    Mutant: mean Mut_Ratio from the two independent mutant lines.
# ---------------------------------------------------------------------------
printf 'panel_id\tgenotype\ttissue\tsite_id\tratio\tn_replicates\n' > "$RATIO_TSV"

for i in "${!tissues[@]}"; do
    tissue="${tissues[$i]}"
    tissue_label="${tissue_labels[$i]}"
    wt_file="$INPUT_DIR/$tissue/$GENOTYPE/WT-${tissue}_sites.Treat_avg_ratio.tsv"
    line1_file="$INPUT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_1-${tissue}.dependent-sites.tsv"
    line2_file="$INPUT_DIR/$tissue/$GENOTYPE/${GENOTYPE}_2-${tissue}.dependent-sites.tsv"
    require_file "$wt_file"
    require_file "$line1_file"
    require_file "$line2_file"

    # WT ratio for this tissue.
    awk -F '\t' -v OFS='\t' \
        -v panel_id="$GENOTYPE" -v tissue_label="$tissue_label" \
        -v out_file="$RATIO_TSV" '
        function read_header(file_number) {
            delete col
            for (j = 1; j <= NF; j++) {
                sub(/\r$/, "", $j)
                col[$j] = j
            }
            if (file_number == 1) {
                required_ok = ("chr" in col) && ("pos" in col) && ("strand" in col) && ("site_id" in col)
            } else {
                required_ok = ("chr" in col) && ("pos" in col) && ("strand" in col) && ("Treat_avg_ratio" in col) && ("Pass" in col)
            }
            if (!required_ok) {
                printf("Required column missing in %s\n", FILENAME) > "/dev/stderr"
                exit 1
            }
        }
        FNR == 1 {
            file_index++
            read_header(file_index)
            next
        }
        file_index == 1 {
            if ($1 == panel_id) union_site[$(col["site_id"])] = 1
            next
        }
        file_index == 2 {
            chr = $(col["chr"]); pos = $(col["pos"]); strand = $(col["strand"])
            input_key = chr "|" pos "|" strand
            if (input_key in wt_ratio) {
                # A duplicate site would silently bias downstream statistics.
                printf("Duplicate site in %s: %s\n", FILENAME, input_key) > "/dev/stderr"
                exit 1
            }
            if ($(col["Pass"]) == 1 && input_key in union_site) {
                wt_ratio[input_key] = $(col["Treat_avg_ratio"]) + 0
            }
            next
        }
        END {
            for (key in wt_ratio) {
                printf "%s\tWT\t%s\t%s\t%.15g\t1\n", panel_id, tissue_label, key, wt_ratio[key] >> out_file
            }
        }
    ' "$UNION_TSV" "$wt_file"

    # Mutant ratio: average the two mutant lines for every union site.
    awk -F '\t' -v OFS='\t' \
        -v panel_id="$GENOTYPE" -v tissue_label="$tissue_label" \
        -v out_file="$RATIO_TSV" '
        function read_header(file_number,    j) {
            delete col
            delete seen_in_file
            for (j = 1; j <= NF; j++) {
                sub(/\r$/, "", $j)
                col[$j] = j
            }
            if (file_number == 1) {
                required_ok = ("chr" in col) && ("pos" in col) && ("strand" in col) && ("site_id" in col)
            } else {
                required_ok = ("chr" in col) && ("pos" in col) && ("strand" in col) && ("Mut_Ratio" in col)
            }
            if (!required_ok) {
                printf("Required column missing in %s\n", FILENAME) > "/dev/stderr"
                exit 1
            }
        }
        function add_ratio(    key) {
            key = $(col["chr"]) "|" $(col["pos"]) "|" $(col["strand"])
            if (key in seen_in_file) {
                printf("Duplicate site within one input file: %s in %s\n", key, FILENAME) > "/dev/stderr"
                exit 1
            }
            seen_in_file[key] = 1
            ratio_sum[key] += $(col["Mut_Ratio"]) + 0
            ratio_n[key]++
        }
        FNR == 1 { file_index++; read_header(file_index); next }
        file_index == 1 {
            if ($1 == panel_id) union_site[$(col["site_id"])] = 1
            next
        }
        file_index > 1 {
            add_ratio()
            next
        }
        END {
            for (key in union_site) {
                if (ratio_n[key] != 2) {
                    printf("Expected two mutant-line ratios for %s, found %d\n", key, ratio_n[key]) > "/dev/stderr"
                    exit 1
                }
                printf "%s\tmutant\t%s\t%s\t%.15g\t%d\n", panel_id, tissue_label, key, ratio_sum[key] / ratio_n[key], ratio_n[key] >> out_file
            }
        }
    ' "$UNION_TSV" "$line1_file" "$line2_file"
done

# Validate that all six groups contain the complete union. This catches silent
# filtering or coordinate-format changes before the table reaches plot.R.
expected_union_n=$(($(wc -l < "$UNION_TSV") - 1))
expected_rows=$((expected_union_n * 6))
actual_rows=$(($(wc -l < "$RATIO_TSV") - 1))
if (( actual_rows != expected_rows )); then
    echo "[data_prep] Incomplete ratio table: expected $expected_rows rows, got $actual_rows." >&2
    exit 1
fi

echo "[data_prep] Generated: $UNION_TSV ($expected_union_n sites)"
echo "[data_prep] Generated: $RATIO_TSV ($actual_rows data rows)"
