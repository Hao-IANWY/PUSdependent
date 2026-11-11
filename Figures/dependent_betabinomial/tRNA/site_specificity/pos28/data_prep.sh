#!/usr/bin/env bash
# data_prep.sh — Sprinzl position 28 site table (Flower panel)
# Question: how many representative (nuclear-encoded) tRNAs are modified at Sprinzl
#           position 28 in WT Flower, and what modification fractions do these tRNAs
#           show in the two alleles each of trua1a and trua5?
# Compartment restricted to Nu (nuclear only; Mt/Pt excluded).
# All logic lives in the shared script ../build_site_specificity_table.sh; this file
# only binds parameters.
# To add a tissue: copy this directory and change TISSUE and the file name; to add a
# position: copy this directory and change POS.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMMON_DIR="$(dirname "$SCRIPT_DIR")"

TISSUE="Flower"
GENOTYPES="trua1a,trua5"
POS=28
COMPARTMENTS="Nu"

bash "$COMMON_DIR/build_site_specificity_table.sh" "$TISSUE" "$GENOTYPES" "$POS" "$SCRIPT_DIR" "$COMPARTMENTS"

test -s "$SCRIPT_DIR/pos${POS}_flower_site_table.tsv"
