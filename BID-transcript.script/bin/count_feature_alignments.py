#!/usr/bin/env python3
"""Count BAM alignments by reference-sequence feature.

Two counting modes are supported:
  * all_alignments: every mapped BAM record contributes one count.
  * primary: only mapped, non-secondary, non-supplementary records contribute.

The output is a raw integer TSV with one row per FASTA feature and one column
per input BAM. Features absent from a BAM retain a zero count.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shlex
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import pysam


VALID_MODES = ("all_alignments", "primary")
_WORKER_FEATURE_IDS: tuple[str, ...] | None = None
_WORKER_MODE: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate a feature-by-sample raw alignment-count matrix."
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Reference FASTA; the first whitespace-delimited header token is the feature ID.",
    )
    parser.add_argument(
        "--mode",
        required=True,
        choices=VALID_MODES,
        help="all_alignments counts every mapped record; primary excludes secondary/supplementary records.",
    )
    parser.add_argument(
        "--bam",
        action="append",
        default=[],
        metavar="SAMPLE=BAM",
        help="Input BAM as sample-name=path. Repeat this option once per sample.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV matrix.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Maximum number of BAMs to count in parallel (default: 1).",
    )
    return parser.parse_args()


def parse_sample_bams(values: list[str]) -> list[tuple[str, str]]:
    parsed: list[tuple[str, str]] = []
    seen: set[str] = set()

    for value in values:
        sample_name, separator, bam_path = value.partition("=")
        if not separator or not sample_name or not bam_path:
            raise ValueError(
                f"Invalid --bam value {shlex.quote(value)!r}; expected SAMPLE=BAM."
            )
        if sample_name in seen:
            raise ValueError(f"Duplicate sample name: {sample_name}")
        seen.add(sample_name)
        parsed.append((sample_name, bam_path))

    if not parsed:
        raise ValueError("At least one --bam SAMPLE=BAM input is required.")
    return parsed


def read_fasta_features(fasta_path: str) -> tuple[str, ...]:
    features: list[str] = []
    seen: set[str] = set()

    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith(">"):
                continue
            feature_id = line[1:].strip().split(None, 1)[0] if line[1:].strip() else ""
            if not feature_id:
                raise ValueError(f"Empty FASTA header at {fasta_path}:{line_number}")
            if feature_id in seen:
                raise ValueError(f"Duplicate FASTA feature ID: {feature_id}")
            features.append(feature_id)
            seen.add(feature_id)

    if not features:
        raise ValueError(f"No FASTA features found in {fasta_path}")
    return tuple(features)


def initialize_worker(feature_ids: tuple[str, ...], mode: str) -> None:
    global _WORKER_FEATURE_IDS, _WORKER_MODE
    _WORKER_FEATURE_IDS = feature_ids
    _WORKER_MODE = mode


def count_bam(bam_path: str) -> dict[str, int]:
    global _WORKER_FEATURE_IDS, _WORKER_MODE
    feature_ids = _WORKER_FEATURE_IDS
    mode = _WORKER_MODE
    if feature_ids is None or mode is None:
        raise RuntimeError("Worker counters were not initialized.")

    counts = {feature_id: 0 for feature_id in feature_ids}
    unknown_references: set[str] = set()

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for alignment in bam.fetch(until_eof=True):
            if alignment.is_unmapped:
                continue
            if mode == "primary" and (
                alignment.is_secondary or alignment.is_supplementary
            ):
                continue

            reference_name = alignment.reference_name
            if reference_name not in counts:
                unknown_references.add(reference_name)
                continue
            counts[reference_name] += 1

    if unknown_references:
        unknown_list = ", ".join(sorted(unknown_references))
        raise ValueError(
            f"{bam_path} contains reference(s) absent from the FASTA: {unknown_list}"
        )

    return counts


def count_all_bams(
    bam_paths: list[str], feature_ids: tuple[str, ...], mode: str, jobs: int
) -> list[dict[str, int]]:
    if jobs <= 0:
        raise ValueError("--jobs must be a positive integer.")

    if jobs == 1:
        global _WORKER_FEATURE_IDS, _WORKER_MODE
        _WORKER_FEATURE_IDS = feature_ids
        _WORKER_MODE = mode
        return [count_bam(path) for path in bam_paths]

    with ProcessPoolExecutor(
        max_workers=jobs,
        initializer=initialize_worker,
        initargs=(feature_ids, mode),
    ) as executor:
        return list(executor.map(count_bam, bam_paths))


def write_matrix(
    output_path: str,
    feature_ids: tuple[str, ...],
    samples: list[str],
    count_rows: list[dict[str, int]],
) -> None:
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary_path = handle.name
        handle.write("feature_id\t" + "\t".join(samples) + "\n")
        for feature_id in feature_ids:
            values = (str(row[feature_id]) for row in count_rows)
            handle.write(feature_id + "\t" + "\t".join(values) + "\n")

    try:
        os.replace(temporary_path, destination)
        # NamedTemporaryFile creates 0600; align with the surrounding 0664 outputs.
        os.chmod(destination, 0o664)
    except BaseException:
        try:
            os.unlink(temporary_path)
        except FileNotFoundError:
            pass
        raise


def main() -> int:
    args = parse_args()
    try:
        samples_and_bams = parse_sample_bams(args.bam)
        samples = [sample for sample, _ in samples_and_bams]
        bam_paths = [path for _, path in samples_and_bams]
        feature_ids = read_fasta_features(args.fasta)
        count_rows = count_all_bams(bam_paths, feature_ids, args.mode, args.jobs)
        write_matrix(args.output, feature_ids, samples, count_rows)
    except (OSError, ValueError, RuntimeError) as error:
        print(f"count_feature_alignments.py: error: {error}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
