#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Adds a Sprinzl_pos column (canonical 0–76) to a site table, mapping native tRNA positions to functional bases.

Background
----------
The reference tRNA set (/home/HL/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa) used for alignment consists of mature tRNA sequences (including 3′-CCA) from the plantRNA database. Each 'pos' represents 1-based coordinates relative to each tRNA. As D-loops, variable regions, and T-loops vary in length, a single native 'pos' maps to different functional positions across different tRNAs. Therefore, positions are linearly recalibrated into fixed canonical slots based on 17 secondary structure regions:

    canonical = region_offset + round_half_up((region_pos - 1) * (L - 1) / (n - 1)) + 1

'n' = actual segment length, 'L' = fixed length (Sprinzl layout, 17 segments, 77 slots total, 0–76; 0 is an extra 5′ base slot used only by His-GTG). Values are rounded half-up (int(x + 0.5)).

Annotations are sourced from the local 'plantRNA_database_unique.tsv' file (220 unique tRNAs, each with 17 region segment definitions). Each tRNA copy uses its own specific annotation, independent of its anticodon family. For sequences not present in this database (e.g., mRNA transcripts, sn/snoRNA), the original row is preserved, and 'Sprinzl_pos' is left blank.

Usage
-----
    python3 add_sprinzl_pos.py --input tRNA/STAT/fisher/mRNA/trua1a_1_Flower_sites.positive.tsv
    python3 add_sprinzl_pos.py --input tRNA/STAT/dependent_betabinomial/ncRNA/Flower/trua1b/trua1b-Flower.dependent.tsv
    python3 add_sprinzl_pos.py -i in.tsv --chr-col chr --pos-col pos -o out.tsv

Output
------
Defaults to <stem>.Sprinzl<ext> in the input directory (e.g., 'trua1b-Flower.dependent.Sprinzl.tsv').
The output contains all original columns with 'Sprinzl_pos' added immediately following the 'pos' column.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

SPRINZL_COLUMN = "Sprinzl_pos"

# Fixed lengths of the 17 structure regions (5'->3'), summing to 77 = canonical 0-76
REGION_NAMES = [
    "p1", "Accstem_5p", "p8-9", "Dstem_5p", "Dloop", "Dstem_3p",
    "p26", "Acstem_5p", "Acloop", "Acstem_3p", "Vregion",
    "Tstem_5p", "Tloop", "Tstem_3p", "Accstem_3p", "p73", "CCA",
]
FIXED_LENGTHS = [1, 7, 2, 4, 8, 4, 1, 5, 7, 5, 5, 5, 7, 5, 7, 1, 3]
REGION_OFFSETS = []
_offset = 0
for _length in FIXED_LENGTHS:
    REGION_OFFSETS.append(_offset)
    _offset += _length


def log(message: str) -> None:
    print(f"[add_sprinzl] {message}", file=sys.stderr)


def read_tsv(path: str) -> tuple[list[str], list[list[str]]]:
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            raise SystemExit(f"ERROR: empty file: {path}")
        rows = [row for row in reader if row]
    header = [name.lstrip("﻿") for name in header]
    return header, rows


def load_annotation(path: str) -> dict[str, tuple[list[int], int]]:
    """Read the annotation TSV; return {ID: (17 region lengths, total sequence length)}."""
    header, rows = read_tsv(path)
    missing = [name for name in ["ID", "Whole_sequence", *REGION_NAMES] if name not in header]
    if missing:
        raise SystemExit(f"ERROR: annotation {path} is missing columns: {', '.join(missing)}")

    id_idx = header.index("ID")
    seq_idx = header.index("Whole_sequence")
    region_idx = [header.index(name) for name in REGION_NAMES]

    annotation: dict[str, tuple[list[int], int]] = {}
    for line_no, row in enumerate(rows, start=2):
        trna_id = row[id_idx].strip()
        region_lengths = [len(row[idx].strip()) for idx in region_idx]
        whole_length = len(row[seq_idx].strip())
        if sum(region_lengths) != whole_length:
            raise SystemExit(
                f"ERROR: annotation row {line_no} ({trna_id}): "
                f"region lengths sum to {sum(region_lengths)} but Whole_sequence is {whole_length}"
            )
        annotation[trna_id] = (region_lengths, whole_length)
    log(f"annotation: {len(annotation)} tRNAs from {path}")
    return annotation


def to_sprinzl(region_lengths: list[int], pos: int) -> int | None:
    """Map a native pos to canonical 0-76; return None if pos is out of range."""
    start = 1
    for idx, length in enumerate(region_lengths):
        if length == 0:
            continue
        end = start + length - 1
        if pos <= end:
            region_pos = pos - start + 1
            fixed = FIXED_LENGTHS[idx]
            if fixed == 1 or length == 1:
                relative = 1
            else:
                relative = int((region_pos - 1) * (fixed - 1) / (length - 1) + 0.5) + 1
            return REGION_OFFSETS[idx] + relative - 1
        start = end + 1
    return None


def parse_pos(value: str) -> int | None:
    text = value.strip()
    try:
        return int(text)
    except ValueError:
        try:
            number = float(text)
        except ValueError:
            return None
        return int(number) if number.is_integer() else None


def default_output_path(input_path: str) -> str:
    root, ext = os.path.splitext(input_path)
    return f"{root}.Sprinzl{ext or '.tsv'}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Append a Sprinzl_pos column (canonical 0-76) after the pos column of a TSV site table.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Input TSV site table")
    parser.add_argument("--chr-col", default="chr", help="Column holding the tRNA / reference sequence ID")
    parser.add_argument("--pos-col", default="pos", help="Column holding the native pos")
    parser.add_argument(
        "--annotation",
        default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "plantRNA_database_unique.tsv"),
        help="Annotation TSV (220 unique tRNAs)",
    )
    parser.add_argument("-o", "--output", default=None, help="Output path (default <stem>.Sprinzl<ext>)")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with an error on rows with an invalid or out-of-range pos (default: warn and leave the cell blank)",
    )
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        raise SystemExit(f"ERROR: input not found: {args.input}")
    if not os.path.isfile(args.annotation):
        raise SystemExit(f"ERROR: annotation not found: {args.annotation}")

    annotation = load_annotation(args.annotation)
    header, rows = read_tsv(args.input)

    for column in (args.chr_col, args.pos_col):
        if column not in header:
            raise SystemExit(
                f"ERROR: column '{column}' not found in {args.input}; available: {', '.join(header)}"
            )
    if SPRINZL_COLUMN in header:
        raise SystemExit(f"ERROR: column '{SPRINZL_COLUMN}' already exists in {args.input}")

    chr_idx = header.index(args.chr_col)
    pos_idx = header.index(args.pos_col)

    output_path = args.output or default_output_path(args.input)
    output_dir = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_dir, exist_ok=True)

    mapped = 0
    skipped_unknown = 0
    skipped_bad_pos = 0
    unknown_ids: dict[str, int] = {}
    bad_examples: list[str] = []

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header[: pos_idx + 1] + [SPRINZL_COLUMN] + header[pos_idx + 1 :])

        for line_no, row in enumerate(rows, start=2):
            if len(row) < len(header):
                row = row + [""] * (len(header) - len(row))
            trna_id = row[chr_idx].strip()
            info = annotation.get(trna_id)
            sprinzl = ""

            if info is None:
                skipped_unknown += 1
                unknown_ids[trna_id] = unknown_ids.get(trna_id, 0) + 1
            else:
                region_lengths, whole_length = info
                pos = parse_pos(row[pos_idx])
                if pos is None or pos < 1 or pos > whole_length:
                    skipped_bad_pos += 1
                    if len(bad_examples) < 5:
                        bad_examples.append(f"line {line_no}: {trna_id} pos={row[pos_idx]!r}")
                    if args.strict:
                        raise SystemExit(f"ERROR: invalid pos at line {line_no}: {trna_id} {row[pos_idx]!r}")
                else:
                    sprinzl = to_sprinzl(region_lengths, pos)
                    if sprinzl is None:
                        skipped_bad_pos += 1
                        if args.strict:
                            raise SystemExit(f"ERROR: unmapped pos at line {line_no}: {trna_id} {pos}")
                    else:
                        mapped += 1

            writer.writerow(row[: pos_idx + 1] + [sprinzl] + row[pos_idx + 1 :])

    log(f"rows: {len(rows)} | mapped: {mapped} | non-tRNA/unknown: {skipped_unknown} | invalid pos: {skipped_bad_pos}")
    if unknown_ids:
        preview = ", ".join(sorted(unknown_ids)[:5])
        log(f"unknown IDs ({len(unknown_ids)} distinct), e.g. {preview}")
    for example in bad_examples:
        log(f"WARNING invalid pos: {example}")
    log(f"output: {output_path}")
    return 1 if (args.strict and skipped_bad_pos) else 0


if __name__ == "__main__":
    sys.exit(main())
