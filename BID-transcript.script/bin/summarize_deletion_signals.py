#!/usr/bin/env python3
"""Summarize raw cpup gap signals and deletion-bearing reads.

The input is one JSON object per ``--dataset`` argument.  Each object must
contain sample, replicate, file_index, ref_type, mpileup, bam, and log paths.
The output is written atomically as a TSV with one row per file_index/ref_type.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Iterable

OUTPUT_COLUMNS = (
    "sample",
    "replicate",
    "file_index",
    "ref_type",
    "final_reads",
    "deletion_signal_reads",
    "deletion_signal_read_ratio",
    "assessed_bases",
    "gap_signal_bases",
    "gap_signal_total",
    "gap_signal_base_ratio",
)
REQUIRED_DATASET_KEYS = frozenset(OUTPUT_COLUMNS[:4]) | {
    "mpileup",
    "bam",
    "log",
}
REQUIRED_LOG_DATASET_KEYS = frozenset(OUTPUT_COLUMNS[:4]) | {"log"}


class DatasetError(ValueError):
    """Raised for an invalid dataset specification."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--dataset",
        action="append",
        default=[],
        metavar="JSON",
        help=(
            "JSON object with sample, replicate, file_index, ref_type, "
            "mpileup, bam, and log paths; repeat once per file_index/ref_type."
        ),
    )
    parser.add_argument(
        "--log-dataset",
        action="append",
        default=[],
        metavar="JSON",
        help=(
            "JSON object with sample, replicate, file_index, ref_type, and "
            "log path; metrics are read from an already-augmented log TSV."
        ),
    )
    parser.add_argument("--output", required=True, help="Output summary TSV path")
    parser.add_argument(
        "--samtools",
        default="samtools",
        help="samtools executable to use for BAM streaming (default: samtools)",
    )
    parser.add_argument(
        "--samtools-threads",
        type=int,
        default=1,
        help="Threads requested for the summary rule; samtools uses one fewer (default: 1)",
    )
    return parser.parse_args()


def parse_datasets(values: Iterable[str]) -> list[dict[str, str]]:
    datasets: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for value in values:
        try:
            record = json.loads(value)
        except json.JSONDecodeError as error:
            raise DatasetError(f"Invalid dataset JSON {value!r}: {error}") from error

        if not isinstance(record, dict):
            raise DatasetError(f"Dataset JSON must be an object: {value!r}")

        missing = REQUIRED_DATASET_KEYS - set(record)
        extra = set(record) - REQUIRED_DATASET_KEYS
        if missing or extra:
            details = []
            if missing:
                details.append("missing keys: " + ", ".join(sorted(missing)))
            if extra:
                details.append("unknown keys: " + ", ".join(sorted(extra)))
            raise DatasetError(
                f"Invalid dataset JSON {value!r}; " + "; ".join(details)
            )

        parsed = {key: str(record[key]) for key in REQUIRED_DATASET_KEYS}
        if any(not value for value in parsed.values()):
            raise DatasetError(f"Dataset JSON contains an empty value: {value!r}")

        key = (parsed["file_index"], parsed["ref_type"])
        if key in seen:
            raise DatasetError(f"Duplicate file_index/ref_type dataset: {'/'.join(key)}")
        seen.add(key)

        for field in ("mpileup", "bam", "log"):
            if not Path(parsed[field]).is_file():
                raise FileNotFoundError(
                    f"Missing {field} for {'/'.join(key)}: {parsed[field]}"
                )

        datasets.append(parsed)

    if not datasets:
        raise DatasetError("At least one --dataset JSON object is required.")

    return sort_datasets(datasets)


def parse_log_datasets(values: list[str]) -> list[dict[str, str]]:
    datasets: list[dict[str, str]] = []
    seen: set[tuple[str, str]] = set()

    for value in values:
        try:
            record = json.loads(value)
        except json.JSONDecodeError as error:
            raise DatasetError(f"Invalid log-dataset JSON {value!r}: {error}") from error

        if not isinstance(record, dict):
            raise DatasetError(f"Log-dataset JSON must be an object: {value!r}")

        missing = REQUIRED_LOG_DATASET_KEYS - set(record)
        extra = set(record) - REQUIRED_LOG_DATASET_KEYS
        if missing or extra:
            details = []
            if missing:
                details.append("missing keys: " + ", ".join(sorted(missing)))
            if extra:
                details.append("unknown keys: " + ", ".join(sorted(extra)))
            raise DatasetError(
                f"Invalid log-dataset JSON {value!r}; " + "; ".join(details)
            )

        parsed = {key: str(record[key]) for key in REQUIRED_LOG_DATASET_KEYS}
        if any(not value for value in parsed.values()):
            raise DatasetError(f"Log-dataset JSON contains an empty value: {value!r}")

        key = (parsed["file_index"], parsed["ref_type"])
        if key in seen:
            raise DatasetError(f"Duplicate file_index/ref_type dataset: {'/'.join(key)}")
        seen.add(key)

        if not Path(parsed["log"]).is_file():
            raise FileNotFoundError(
                f"Missing log for {'/'.join(key)}: {parsed['log']}"
            )
        datasets.append(parsed)

    if not datasets:
        raise DatasetError("At least one --log-dataset JSON object is required.")
    return sort_datasets(datasets)


def sort_datasets(datasets: list[dict[str, str]]) -> list[dict[str, str]]:
    return sorted(
        datasets,
        key=lambda item: (
            item["sample"],
            item["replicate"],
            item["ref_type"],
            item["file_index"],
        ),
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


def read_final_reads(log_path: str, ref_type: str) -> int:
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

    field = f"{ref_type}_filtered_reads"
    try:
        raw_value = values[headers.index(field)]
    except ValueError as error:
        raise ValueError(f"{log_path}: missing required column {field}.") from error

    try:
        final_reads = int(raw_value)
    except ValueError as error:
        raise ValueError(f"{log_path}: {field} is not an integer: {raw_value!r}.") from error

    if final_reads < 0:
        raise ValueError(f"{log_path}: {field} must be non-negative.")
    return final_reads


def run_pipeline(commands: list[list[str]], label: str) -> str:
    """Run a stdout->stdin command pipeline and return the final stdout.

    Every stage's exit code is checked so an upstream failure (e.g. samtools
    erroring out mid-stream) cannot masquerade as valid downstream output.
    """
    processes = []
    upstream = None
    for argv in commands:
        process = subprocess.Popen(
            argv,
            stdin=upstream,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        if upstream is not None:
            upstream.close()
        upstream = process.stdout
        processes.append(process)

    stdout, last_stderr = processes[-1].communicate()
    errors = []
    for argv, process in zip(commands[:-1], processes[:-1]):
        stderr = process.stderr.read() if process.stderr is not None else ""
        return_code = process.wait()
        if return_code != 0:
            errors.append(f"{' '.join(argv)} exited {return_code}: {stderr.strip()}")
    if processes[-1].returncode != 0:
        errors.append(
            f"{' '.join(commands[-1])} exited {processes[-1].returncode}: "
            f"{last_stderr.strip()}"
        )
    if errors:
        raise RuntimeError(f"{label} failed: " + "; ".join(errors))
    return stdout.strip()


# cpup TSV 是程序生成的固定格式，这里只保留字段数校验；
# 用 awk 单遍聚合代替 Python 逐行解析，GB 级文件快一个量级。
MPILEUP_AWK_PROGRAM = r"""
BEGIN { FS = "\t" }
NF != 5 {
    printf "expected 5 tab-separated fields, found %d (line %d)\n", NF, NR > "/dev/stderr"
    exit 1
}
$4 != "+" { next }
{
    n = split($5, counts, ",")
    if (n != 12) {
        printf "expected 12 cpup count fields, found %d (line %d)\n", n, NR > "/dev/stderr"
        exit 1
    }
    # cpup 的 depth 列只计非 gap 碱基，加上 gap 占位才是该位点的总观测深度。
    observed = counts[1] + counts[8]
    if (observed > 0) {
        assessed++
        if (counts[8] > 0) {
            gap_bases++
            gap_total += counts[8]
        }
    }
}
END { printf "%d\t%d\t%d\n", assessed + 0, gap_bases + 0, gap_total + 0 }
"""


def summarize_mpileup(path: str) -> tuple[int, int, int]:
    """Return covered bases, gap-bearing covered bases, and total gap signals."""
    if not Path(path).is_file():
        raise FileNotFoundError(f"mpileup TSV not found: {path}")
    output = run_pipeline(
        [["awk", MPILEUP_AWK_PROGRAM, path]],
        f"awk mpileup-summary {path}",
    )
    fields = output.split("\t")
    if len(fields) != 3:
        raise ValueError(f"awk mpileup-summary returned unexpected output: {output!r}")
    return int(fields[0]), int(fields[1]), int(fields[2])


def count_deletion_reads(
    bam_path: str,
    samtools: str = "samtools",
    threads: int = 1,
) -> int:
    if threads < 1:
        raise ValueError("--samtools-threads must be at least 1.")
    if not Path(bam_path).is_file():
        raise FileNotFoundError(f"BAM not found: {bam_path}")

    # samtools 负责 BAM 解码，awk 统计 CIGAR 含 D 的唯一 QNAME；
    # 比 Python 逐行 split 快一个量级，结果完全一致。
    output = run_pipeline(
        [
            [samtools, "view", "-@", str(threads - 1), bam_path],
            ["awk", '$6 ~ /D/ && !seen[$1]++ { count++ } END { print count + 0 }'],
        ],
        f"samtools view {bam_path} | awk deletion-count",
    )
    return int(output)


def format_ratio(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "NA"
    return f"{numerator / denominator:.6f}"


def summarize_log_dataset(dataset: dict[str, str]) -> dict[str, int | str]:
    row = read_log_row(dataset["log"])
    ref_type = dataset["ref_type"]
    fields = {
        "final_reads": f"{ref_type}_filtered_reads",
        "deletion_signal_reads": f"{ref_type}_deletion_signal_reads",
        "deletion_signal_read_ratio": f"{ref_type}_deletion_signal_read_ratio",
        "assessed_bases": f"{ref_type}_assessed_bases",
        "gap_signal_bases": f"{ref_type}_gap_signal_bases",
        "gap_signal_total": f"{ref_type}_gap_signal_total",
        "gap_signal_base_ratio": f"{ref_type}_gap_signal_base_ratio",
    }
    missing = [column for column in fields.values() if column not in row]
    if missing:
        raise ValueError(
            f"{dataset['log']}: missing deletion column(s): " + ", ".join(missing)
        )

    parsed: dict[str, int | str] = {
        "sample": dataset["sample"],
        "replicate": dataset["replicate"],
        "file_index": dataset["file_index"],
        "ref_type": ref_type,
    }
    integer_fields = (
        "final_reads",
        "deletion_signal_reads",
        "assessed_bases",
        "gap_signal_bases",
        "gap_signal_total",
    )
    for output_name in integer_fields:
        column = fields[output_name]
        try:
            value = int(row[column])
        except ValueError as error:
            raise ValueError(
                f"{dataset['log']}: {column} is not an integer: {row[column]!r}."
            ) from error
        if value < 0:
            raise ValueError(f"{dataset['log']}: {column} must be non-negative.")
        parsed[output_name] = value

    for output_name in (
        "deletion_signal_read_ratio",
        "gap_signal_base_ratio",
    ):
        column = fields[output_name]
        value = row[column]
        if value != "NA":
            try:
                ratio = float(value)
            except ValueError as error:
                raise ValueError(
                    f"{dataset['log']}: {column} is not numeric or NA: {value!r}."
                ) from error
            if not 0.0 <= ratio <= 1.0:
                raise ValueError(f"{dataset['log']}: {column} is outside [0, 1].")
        parsed[output_name] = value

    if int(parsed["deletion_signal_reads"]) > int(parsed["final_reads"]):
        raise ValueError(
            f"{dataset['log']}: deletion-bearing reads exceed filtered reads."
        )
    if int(parsed["gap_signal_bases"]) > int(parsed["assessed_bases"]):
        raise ValueError(
            f"{dataset['log']}: gap-bearing bases exceed assessed bases."
        )
    return parsed


def summarize_dataset(
    dataset: dict[str, str], samtools: str, threads: int
) -> dict[str, int | str]:
    final_reads = read_final_reads(dataset["log"], dataset["ref_type"])
    deletion_reads = count_deletion_reads(
        dataset["bam"], samtools=samtools, threads=threads
    )
    if deletion_reads > final_reads:
        raise ValueError(
            f"{dataset['bam']}: unique deletion-bearing QNAMEs ({deletion_reads}) "
            f"exceed filtered alignments ({final_reads})."
        )

    assessed_bases, gap_bases, gap_total = summarize_mpileup(dataset["mpileup"])
    return {
        "sample": dataset["sample"],
        "replicate": dataset["replicate"],
        "file_index": dataset["file_index"],
        "ref_type": dataset["ref_type"],
        "final_reads": final_reads,
        "deletion_signal_reads": deletion_reads,
        "deletion_signal_read_ratio": format_ratio(deletion_reads, final_reads),
        "assessed_bases": assessed_bases,
        "gap_signal_bases": gap_bases,
        "gap_signal_total": gap_total,
        "gap_signal_base_ratio": format_ratio(gap_bases, assessed_bases),
    }


def write_summary(rows: Iterable[dict[str, int | str]], output_path: str) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")

    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
            for row in rows:
                handle.write("\t".join(str(row[column]) for column in OUTPUT_COLUMNS) + "\n")
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    if bool(args.dataset) == bool(args.log_dataset):
        raise DatasetError(
            "Provide either --dataset records or --log-dataset records, not both."
        )

    rows = []
    if args.dataset:
        datasets = parse_datasets(args.dataset)
        for index, dataset in enumerate(datasets, start=1):
            label = f"{dataset['file_index']}/{dataset['ref_type']}"
            print(f"[deletion-stats] {index}/{len(datasets)} {label}", file=sys.stderr)
            rows.append(
                summarize_dataset(dataset, args.samtools, args.samtools_threads)
            )
    else:
        datasets = parse_log_datasets(args.log_dataset)
        for index, dataset in enumerate(datasets, start=1):
            label = f"{dataset['file_index']}/{dataset['ref_type']}"
            print(f"[deletion-stats] {index}/{len(datasets)} {label} from log", file=sys.stderr)
            rows.append(summarize_log_dataset(dataset))

    write_summary(rows, args.output)
    print(f"[deletion-stats] Wrote {len(rows)} rows to {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (DatasetError, OSError, ValueError, subprocess.CalledProcessError) as error:
        print(f"[deletion-stats] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
