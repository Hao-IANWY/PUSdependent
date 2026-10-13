#!/usr/bin/env python3
"""Collate per-library mapping logs into a long processing-summary TSV."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Iterable

OUTPUT_COLUMNS = ("Sample", "Replicate", "Stage", "Reads")
REQUIRED_KEYS = frozenset(("sample", "replicate", "file_index", "log"))


class DatasetError(ValueError):
    """Raised for an invalid dataset specification."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        metavar="JSON",
        help=(
            "JSON object with sample, replicate, file_index, and log; "
            "repeat once per library."
        ),
    )
    parser.add_argument(
        "--ref-types",
        required=True,
        help="Comma-separated reference types in processing-stage order",
    )
    parser.add_argument("--output", required=True, help="Output long-format TSV")
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    fields = [field.strip() for field in value.split(",")]
    if not fields or any(not field for field in fields):
        raise ValueError("--ref-types must contain non-empty comma-separated names.")
    if len(set(fields)) != len(fields):
        raise ValueError("--ref-types contains duplicate names.")
    return fields


def parse_datasets(values: list[str]) -> list[dict[str, str]]:
    datasets: list[dict[str, str]] = []
    seen: set[str] = set()

    for value in values:
        try:
            record = json.loads(value)
        except json.JSONDecodeError as error:
            raise DatasetError(f"Invalid dataset JSON {value!r}: {error}") from error
        if not isinstance(record, dict):
            raise DatasetError(f"Dataset JSON must be an object: {value!r}")

        missing = REQUIRED_KEYS - set(record)
        extra = set(record) - REQUIRED_KEYS
        if missing or extra:
            details = []
            if missing:
                details.append("missing keys: " + ", ".join(sorted(missing)))
            if extra:
                details.append("unknown keys: " + ", ".join(sorted(extra)))
            raise DatasetError(f"Invalid dataset JSON {value!r}; " + "; ".join(details))

        parsed = {key: str(record[key]) for key in REQUIRED_KEYS}
        if any(not value for value in parsed.values()):
            raise DatasetError(f"Dataset JSON contains an empty value: {value!r}")
        if parsed["file_index"] in seen:
            raise DatasetError(f"Duplicate file_index: {parsed['file_index']}")
        seen.add(parsed["file_index"])
        if not Path(parsed["log"]).is_file():
            raise FileNotFoundError(
                f"Missing log for {parsed['file_index']}: {parsed['log']}"
            )
        datasets.append(parsed)

    if not datasets:
        raise DatasetError("At least one --dataset JSON object is required.")
    return sorted(
        datasets,
        key=lambda item: (item["sample"], item["replicate"], item["file_index"]),
    )


def read_log_row(log_path: str) -> dict[str, str]:
    with open(log_path, "r", encoding="utf-8", newline="") as handle:
        header_line = handle.readline()
        value_line = handle.readline()
        if handle.readline().strip():
            raise ValueError(f"{log_path}: expected exactly one data row.")
    if not header_line or not value_line:
        raise ValueError(f"{log_path}: missing header or data row.")

    headers = header_line.rstrip("\r\n").split("\t")
    values = value_line.rstrip("\r\n").split("\t")
    if len(headers) != len(values):
        raise ValueError(
            f"{log_path}: header has {len(headers)} fields, data has {len(values)}."
        )
    return dict(zip(headers, values))


def parse_read_count(log_path: str, field: str, row: dict[str, str]) -> int:
    if field not in row:
        raise ValueError(f"{log_path}: missing processing column {field}.")
    try:
        value = int(row[field])
    except ValueError as error:
        raise ValueError(f"{log_path}: {field} is not an integer: {row[field]!r}.") from error
    if value < 0:
        raise ValueError(f"{log_path}: {field} must be non-negative.")
    return value


def build_rows(
    datasets: Iterable[dict[str, str]], ref_types: list[str]
) -> list[tuple[str, str, str, int]]:
    rows = []
    fixed_stages = (
        ("Raw_Reads", "raw_reads"),
        ("Trimmed_Reads", "trimmed_reads"),
        ("Dedup_Reads", "dedup_reads"),
    )
    for dataset in datasets:
        row = read_log_row(dataset["log"])
        for stage, field in fixed_stages:
            rows.append(
                (dataset["sample"], dataset["replicate"], stage, parse_read_count(dataset["log"], field, row))
            )
        for ref_type in ref_types:
            for suffix in ("unique", "multi"):
                field = f"{ref_type}_{suffix}"
                rows.append(
                    (
                        dataset["sample"],
                        dataset["replicate"],
                        f"{ref_type}_{suffix}",
                        parse_read_count(dataset["log"], field, row),
                    )
                )
    return rows


def write_rows(rows: list[tuple[str, str, str, int]], output_path: str) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
            for row in rows:
                handle.write("\t".join(map(str, row)) + "\n")
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    ref_types = split_csv(args.ref_types)
    datasets = parse_datasets(args.dataset)
    rows = build_rows(datasets, ref_types)
    write_rows(rows, args.output)
    print(f"[processing-stats] Wrote {len(rows)} rows to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (DatasetError, OSError, ValueError) as error:
        print(f"[processing-stats] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
