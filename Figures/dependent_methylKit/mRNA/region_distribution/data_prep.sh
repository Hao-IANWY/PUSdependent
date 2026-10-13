#!/usr/bin/env bash
# data_prep.sh — Map mRNA-dependent transcript sites to genome coordinates and
# annotate their transcript regions with Araport11_fixed.gtf.
#
# Inputs:
#   polyARNA/STAT/dependent_methylKit/mRNA/<tissue>/<genotype>/<genotype>-<tissue>.dependent.tsv
#   Araport11_fixed.gtf
# Outputs:
#   site_analyse/<genotype>-<tissue>.dependent.{bed,annoated.bed,annoated.tsv}

set -euo pipefail
IFS=$'\n\t'

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT_DIR="$(cd "$SCRIPT_DIR/../../../.." && pwd)"
readonly INPUT_DIR="$ROOT_DIR/polyARNA/STAT/dependent_methylKit/mRNA"
readonly OUTPUT_DIR="$SCRIPT_DIR/site_analyse"
readonly PYTHON_BIN="/home/HL/miniforge3/envs/analysis/bin/python3.11"
readonly GPPY="/home/HL/software/gppy/gppy/gtf.py"
readonly GTF="$ROOT_DIR/Araport11_fixed.gtf"
readonly ANNOTATION_SCRIPT="/home/haoli/software/Stript/BASAL_PARISE/bed_annoat.R"

readonly -a TISSUES=(Root leaf Flower)
readonly -a GENOTYPES=(trua1a trua1b trua2 trua3 trua4 trua5 trua6)

log() {
    printf '[data_prep] %s\n' "$*" >&2
}

error() {
    log "ERROR: $*"
    exit 1
}

require_command() {
    local cmd="$1"
    command -v "$cmd" >/dev/null 2>&1 || error "Required command not found: $cmd"
}

require_nonempty_file() {
    local path="$1"
    [[ -s "$path" ]] || error "Required file not found or empty: $path"
}

line_count() {
    local path="$1"
    if [[ -s "$path" ]]; then
        wc -l < "$path" | awk '{print $1 + 0}'
    else
        printf '0\n'
    fi
}

annotate_one_file() {
    local tissue="$1"
    local genotype="$2"
    local input_file="$INPUT_DIR/$tissue/$genotype/$genotype-$tissue.dependent.tsv"
    local base_name="$genotype-$tissue.dependent"
    local body_file="$OUTPUT_DIR/$base_name.tsv"
    local bed_file="$OUTPUT_DIR/$base_name.bed"
    local annotated_bed_file="$OUTPUT_DIR/$base_name.annoated.bed"
    local annotated_tsv_file="$OUTPUT_DIR/$base_name.annoated.tsv"

    require_nonempty_file "$input_file"  # header-only files are valid zero-site results
    log "Processing $tissue / $genotype"

    # Normalize CRLF endings while retaining the original column content.
    awk 'NR > 1 { sub(/\r$/, ""); print }' "$input_file" > "$body_file"

    # methylKit can call no sites for a genotype/tissue. Preserve the schema so
    # downstream code still sees all 21 expected panels.
    if [[ ! -s "$body_file" ]]; then
        awk 'NR == 1 { sub(/\r$/, ""); print $0 "\tannoation" }' \
            "$input_file" > "$annotated_tsv_file"
        log "Wrote $annotated_tsv_file (0 sites)"
        return 0
    fi
    local input_count body_count
    input_count=$(( $(wc -l < "$input_file") - 1 ))
    body_count=$(line_count "$body_file")
    [[ "$input_count" -eq "$body_count" ]] || \
        error "Row-count mismatch for $input_file: expected $input_count, wrote $body_count"

    # gtf.py t2g maps transcript coordinates to Araport11 genome coordinates.
    # Its last three output columns are chromosome, strand, and genomic position.
    "$PYTHON_BIN" "$GPPY" t2g -g "$GTF" -i "$body_file" |
        awk 'BEGIN { FS = "\t"; OFS = "\t"; bad = 0 }
        {
            if ($(NF - 2) == "NA" || $(NF - 1) == "NA" || $NF == "NA") {
                printf("No GTF mapping for %s:%s:%s\n", $1, $2, $3) > "/dev/stderr"
                bad = 1
                next
            }
            print $(NF - 2), $NF - 1, $NF, ".", ".", $(NF - 1)
        }
        END { if (bad) exit 1 }' > "$bed_file"
    require_nonempty_file "$bed_file"

    local bed_count
    bed_count=$(line_count "$bed_file")
    [[ "$body_count" -eq "$bed_count" ]] || \
        error "BED row-count mismatch for $input_file: expected $body_count, found $bed_count"

    Rscript "$ANNOTATION_SCRIPT" \
        -b "$bed_file" \
        -g "$GTF" \
        -o "$annotated_bed_file" >&2
    require_nonempty_file "$annotated_bed_file"

    local annotated_bed_count
    annotated_bed_count=$(line_count "$annotated_bed_file")
    [[ "$body_count" -eq "$annotated_bed_count" ]] || \
        error "Annotation row-count mismatch for $input_file: expected $body_count, found $annotated_bed_count"

    {
        awk 'NR == 1 { sub(/\r$/, ""); print $0 "\tannoation" }' "$input_file"
        paste -d '\t' "$body_file" \
            <(awk -F '\t' '{ sub(/\r$/, ""); print $4 }' "$annotated_bed_file")
    } > "$annotated_tsv_file"
    require_nonempty_file "$annotated_tsv_file"

    local annotated_count
    annotated_count=$(( $(wc -l < "$annotated_tsv_file") - 1 ))
    [[ "$body_count" -eq "$annotated_count" ]] || \
        error "Annotated TSV row-count mismatch for $input_file: expected $body_count, found $annotated_count"

    log "Wrote $annotated_tsv_file ($annotated_count sites)"
}

main() {
    require_command "$PYTHON_BIN"
    require_command Rscript
    require_command awk
    require_command paste

    require_nonempty_file "$GPPY"
    require_nonempty_file "$GTF"
    require_nonempty_file "$ANNOTATION_SCRIPT"

    mkdir -p "$OUTPUT_DIR"
    find "$OUTPUT_DIR" -maxdepth 1 -type f -delete

    local tissue genotype annotated_count expected_count
    annotated_count=0
    expected_count=0
    for tissue in "${TISSUES[@]}"; do
        for genotype in "${GENOTYPES[@]}"; do
            annotate_one_file "$tissue" "$genotype"
            expected_count=$((expected_count + 1))
            if [[ -s "$OUTPUT_DIR/$genotype-$tissue.dependent.annoated.tsv" ]]; then
                annotated_count=$((annotated_count + 1))
            fi
        done
    done

    [[ "$annotated_count" -eq "$expected_count" ]] || \
        error "Expected $expected_count annotated files, found $annotated_count"

    log "All $annotated_count mRNA dependent-region annotations completed."
}

main "$@"
