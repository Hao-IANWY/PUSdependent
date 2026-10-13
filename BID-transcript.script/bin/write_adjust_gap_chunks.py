#!/usr/bin/env python3
"""Build deterministic reference-sequence chunks for parallel adjustGap."""

import argparse
import os
from pathlib import Path
from typing import List, Tuple


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fai", required=True, help="Reference FASTA index")
    parser.add_argument("--output", required=True, help="Output chunk manifest TSV")
    parser.add_argument("--chunk-bp", type=int, required=True,
                        help="Close a chunk after at least this many reference bases")
    parser.add_argument("--max-sequences", type=int, required=True,
                        help="Maximum number of reference sequences per chunk")
    return parser.parse_args()


def read_fai(path: Path) -> List[Tuple[str, int]]:
    records: List[Tuple[str, int]] = []
    seen = set()
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1].isdigit():
                raise ValueError(f"{path}:{line_number}: invalid FAI record: {line.rstrip()!r}")
            name, length = fields[0], int(fields[1])
            if name in seen:
                raise ValueError(f"{path}:{line_number}: duplicate sequence name: {name}")
            if length <= 0:
                raise ValueError(f"{path}:{line_number}: non-positive sequence length: {name}")
            seen.add(name)
            records.append((name, length))
    if not records:
        raise ValueError(f"{path}: no reference records")
    return records


def build_chunks(records: List[Tuple[str, int]], chunk_bp: int, max_sequences: int) -> List[List[str]]:
    chunks: List[List[str]] = []
    current: List[str] = []
    current_bp = 0

    def close_chunk() -> None:
        nonlocal current, current_bp
        if current:
            chunks.append(current)
        current = []
        current_bp = 0

    for _name, length in records:
        if current and (current_bp >= chunk_bp or len(current) >= max_sequences):
            close_chunk()
        current.append(_name)
        current_bp += length

    close_chunk()
    return chunks


def main() -> None:
    args = parse_args()
    if args.chunk_bp <= 0:
        raise ValueError("--chunk-bp must be a positive integer")
    if args.max_sequences <= 0:
        raise ValueError("--max-sequences must be a positive integer")

    records = read_fai(Path(args.fai))
    chunks = build_chunks(records, args.chunk_bp, args.max_sequences)
    width = max(6, len(str(len(chunks))))
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_output = output.with_name(f".{output.name}.tmp.{os.getpid()}")

    try:
        sequence_lengths = dict(records)
        with temp_output.open("w", encoding="utf-8", newline="") as handle:
            handle.write("chunk_id\tchrom\tlength\n")
            for index, sequences in enumerate(chunks, start=1):
                chunk_id = f"chunk_{index:0{width}d}"
                for name in sequences:
                    handle.write(f"{chunk_id}\t{name}\t{sequence_lengths[name]}\n")
        os.replace(temp_output, output)
    except BaseException:
        temp_output.unlink(missing_ok=True)
        raise

    print(
        f"Wrote {len(chunks)} adjustGap chunks covering {len(records)} reference sequences "
        f"to {output}"
    )


if __name__ == "__main__":
    main()
