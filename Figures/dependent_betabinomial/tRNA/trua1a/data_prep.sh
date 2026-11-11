#!/usr/bin/env bash
# data_prep.sh — Build this genotype's tRNA writer-dependent Ψ site tables (per-tissue + cross-tissue union)
# Purpose: data only, no plotting. All logic lives in the shared script
# ../build_genotype_dependent_tables.sh; this file only binds the genotype and the
# output directory.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENOTYPE="trua1a"

bash "$SCRIPT_DIR/../build_genotype_dependent_tables.sh" "$GENOTYPE" "$SCRIPT_DIR"

test -s "$SCRIPT_DIR/${GENOTYPE}_all_tissues.dependent.tsv"
test -s "$SCRIPT_DIR/${GENOTYPE}_union.dependent.tsv"
