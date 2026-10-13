#!/usr/bin/env python3
"""Build one per-library mapping/deletion log TSV."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from concurrent.futures import ThreadPoolExecutor
import gzip
import importlib.util
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from typing import Iterator

SCRIPT_DIR = Path(__file__).resolve().parent
DELETION_MODULE_PATH = SCRIPT_DIR / "summarize_deletion_signals.py"


def load_deletion_module():
    spec = importlib.util.spec_from_file_location(
        "summarize_deletion_signals", DELETION_MODULE_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {DELETION_MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


deletion = load_deletion_module()


def positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--samtools-threads must be an integer.") from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("--samtools-threads must be at least 1.")
    return parsed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--file-index", required=True)
    parser.add_argument("--read-layout", required=True, choices=("paired", "single"))
    parser.add_argument("--rmdup", required=True, choices=("true", "false"))
    parser.add_argument("--raw-fq", required=True)
    parser.add_argument("--trimmed-fq", required=True)
    parser.add_argument("--dedup-fq", required=True)
    parser.add_argument("--bowtie2-log", required=True)
    parser.add_argument("--star-log")
    parser.add_argument(
        "--bam",
        action="append",
        default=[],
        metavar="REF_TYPE,STAGE,BAM",
        help="Reference type, BAM stage, and path; stage is forward/realign/filtered.",
    )
    parser.add_argument(
        "--mpileup",
        action="append",
        default=[],
        metavar="REF_TYPE,TSV",
        help="Reference type and raw cpup mpileup TSV path.",
    )
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument(
        "--samtools-threads",
        type=positive_int,
        default=1,
        help=(
            "Thread budget for this job; the independent counters run "
            "concurrently with this many workers."
        ),
    )
    return parser.parse_args()


def parse_key_value_args(values: list[str], field_count: int, option: str):
    parsed = []
    seen = set()
    for value in values:
        fields = value.split(",")
        if len(fields) != field_count or any(not field for field in fields):
            raise ValueError(
                f"Invalid {option} value {value!r}; expected {field_count} comma-separated fields."
            )
        key = tuple(fields[:-1])
        if key in seen:
            raise ValueError(f"Duplicate {option} key: {'/'.join(key)}")
        seen.add(key)
        parsed.append((key, fields[-1]))
    return parsed


def read_chunks(path: str) -> Iterator[bytes]:
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            yield chunk


def gzip_chunks(path: str) -> Iterator[bytes]:
    with gzip.open(path, "rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            yield chunk


def count_fastq_reads(path: str) -> int:
    if not Path(path).is_file():
        raise FileNotFoundError(f"FASTQ not found: {path}")
    if str(path).endswith(".gz") and shutil.which("pigz"):
        # pigz 多线程解压比 Python gzip 模块快数倍，wc -l 只负责数行。
        output = deletion.run_pipeline(
            [["pigz", "-dc", str(path)], ["wc", "-l"]],
            f"pigz -dc {path} | wc -l",
        )
        lines = int(output)
    else:
        chunks = gzip_chunks(path) if str(path).endswith(".gz") else read_chunks(path)
        lines = sum(chunk.count(b"\n") for chunk in chunks)
    if lines % 4 != 0:
        raise ValueError(f"FASTQ line count is not divisible by 4: {path}")
    return lines // 4


def run_output(command: list[str], label: str) -> str:
    if Path(command[0]).is_file() or "/" in command[0]:
        executable = command[0]
    else:
        executable = str(Path(command[0]))
    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{label} failed with exit code {completed.returncode}: {completed.stderr.strip()}"
        )
    return completed.stdout.strip()


def count_bam_reads(path: str, samtools: str, threads: int) -> int:
    if not Path(path).is_file():
        raise FileNotFoundError(f"BAM not found: {path}")
    output = run_output(
        [samtools, "view", "-c", "-@", str(max(threads - 1, 0)), path],
        f"samtools view -c {path}",
    )
    try:
        count = int(output)
    except ValueError as error:
        raise ValueError(f"samtools returned non-integer count for {path}: {output!r}") from error
    if count < 0:
        raise ValueError(f"Negative BAM count for {path}: {count}")
    return count


def parse_bowtie2_log(path: str) -> tuple[int, int]:
    unique = 0
    multi = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            # bowtie2 indents summary lines; re.match anchors at position 0.
            line = line.strip()
            match_unique = re.match(r"(\d+)\s+\(.+\)\s+aligned exactly 1 time", line)
            if match_unique:
                unique = int(match_unique.group(1))
                continue
            match_multi = re.match(r"(\d+)\s+\(.+\)\s+aligned >1 times", line)
            if match_multi:
                multi = int(match_multi.group(1))
    return unique, multi


def parse_star_log(path: str) -> tuple[int, int]:
    unique = 0
    multi = 0
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            if "Uniquely mapped reads number" in line:
                unique = int(re.search(r"\d+", line).group())
            elif "Number of reads mapped to multiple loci" in line:
                multi = int(re.search(r"\d+", line).group())
            elif "Number of reads mapped to too many loci" in line:
                multi += int(re.search(r"\d+", line).group())
    return unique, multi


def validate_paths(paths: list[str]) -> None:
    missing = [path for path in paths if not Path(path).is_file()]
    if missing:
        raise FileNotFoundError("Missing input(s): " + ", ".join(missing))


def write_log(output_path: str, stats: "OrderedDict[str, object]") -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.tmp.{os.getpid()}")
    try:
        with temporary.open("w", encoding="utf-8", newline="") as handle:
            handle.write("\t".join(stats.keys()) + "\n")
            handle.write("\t".join(map(str, stats.values())) + "\n")
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def main() -> int:
    args = parse_args()
    bam_records = parse_key_value_args(args.bam, 3, "--bam")
    mpileup_records = parse_key_value_args(args.mpileup, 2, "--mpileup")

    ref_types = []
    for (ref_type, _stage), _path in bam_records:
        if ref_type not in ref_types:
            ref_types.append(ref_type)
    mpileup_ref_types = list(
        dict.fromkeys(key[0] for key, _path in mpileup_records)
    )
    if ref_types != mpileup_ref_types:
        raise ValueError("--bam and --mpileup reference types do not match.")
    if not ref_types:
        raise ValueError("At least one --bam/--mpileup reference type is required.")

    expected_stages = {"forward", "realign", "filtered"}
    stage_paths: dict[tuple[str, str], str] = {}
    for (ref_type, stage), path in bam_records:
        if stage not in expected_stages:
            raise ValueError(f"Invalid BAM stage {stage!r}; expected one of {sorted(expected_stages)}.")
        stage_paths[(ref_type, stage)] = path
    for ref_type in ref_types:
        for stage in expected_stages:
            if (ref_type, stage) not in stage_paths:
                raise ValueError(f"Missing --bam entry for {ref_type}/{stage}.")

    # parse_key_value_args keys are tuples; --mpileup has a single key field, so
    # flatten ("ncRNA",) -> "ncRNA" for per-ref-type lookup below.
    mpileup_paths = {key[0]: path for key, path in mpileup_records}
    validate_paths(
        [
            args.raw_fq,
            args.trimmed_fq,
            args.dedup_fq,
            args.bowtie2_log,
            *(path for path in stage_paths.values()),
            *(path for path in mpileup_paths.values()),
            *( [args.star_log] if args.star_log else [] ),
        ]
    )

    # 所有计数任务彼此独立（3 个 FASTQ、3x ref_type 个 BAM 计数、
    # 2x ref_type 个缺失统计），用线程池并发执行；任务都是 subprocess/IO 型，
    # 单 job 墙钟约等于最慢的一个任务。每个 samtools 调用只带 1 个解压线程，
    # 避免线程池 x 进程内线程的过度超订。
    per_call_samtools_threads = 2
    futures = {}
    with ThreadPoolExecutor(max_workers=args.samtools_threads) as pool:
        futures["raw_reads"] = pool.submit(count_fastq_reads, args.raw_fq)
        futures["trimmed_reads"] = pool.submit(count_fastq_reads, args.trimmed_fq)
        futures["dedup_reads"] = pool.submit(count_fastq_reads, args.dedup_fq)
        for ref_type in ref_types:
            for stage in ("forward", "realign", "filtered"):
                futures[f"{ref_type}_{stage}_reads"] = pool.submit(
                    count_bam_reads,
                    stage_paths[(ref_type, stage)],
                    args.samtools,
                    per_call_samtools_threads,
                )
        for ref_type in ref_types:
            futures[f"{ref_type}_deletion_signal_reads"] = pool.submit(
                deletion.count_deletion_reads,
                stage_paths[(ref_type, "filtered")],
                args.samtools,
                per_call_samtools_threads,
            )
            futures[f"{ref_type}_mpileup_summary"] = pool.submit(
                deletion.summarize_mpileup, mpileup_paths[ref_type]
            )
        counts = {key: future.result() for key, future in futures.items()}

    bowtie_unique, bowtie_multi = parse_bowtie2_log(args.bowtie2_log)
    # Column names follow the actual reference types: Bowtie2 always maps to the
    # first reference, STAR (when configured) maps the unmapped reads to the
    # second. summarize_processing_logs.py expects <ref>_unique/<ref>_multi for
    # every reference type, so the STAR pair is emitted as zeros when the log is
    # unavailable.
    bowtie_ref = ref_types[0]
    star_ref = ref_types[1] if len(ref_types) > 1 else None
    if args.star_log and star_ref is None:
        raise ValueError(
            "--star-log was provided but only one reference type is configured."
        )

    stats: OrderedDict[str, object] = OrderedDict(
        (
            ("read_layout", args.read_layout),
            ("rmdup_enabled", args.rmdup),
            ("raw_reads", counts["raw_reads"]),
            ("trimmed_reads", counts["trimmed_reads"]),
            ("dedup_reads", counts["dedup_reads"]),
            (f"{bowtie_ref}_unique", bowtie_unique),
            (f"{bowtie_ref}_multi", bowtie_multi),
        )
    )
    if star_ref is not None:
        if args.star_log:
            star_unique, star_multi = parse_star_log(args.star_log)
        else:
            star_unique, star_multi = 0, 0
        stats[f"{star_ref}_unique"] = star_unique
        stats[f"{star_ref}_multi"] = star_multi

    for ref_type in ref_types:
        for stage in ("forward", "realign", "filtered"):
            stats[f"{ref_type}_{stage}_reads"] = counts[f"{ref_type}_{stage}_reads"]
            print(
                f"[library-log] {args.file_index} {ref_type} {stage}: "
                f"{stats[f'{ref_type}_{stage}_reads']}",
                file=sys.stderr,
            )

    for ref_type in ref_types:
        filtered_reads = counts[f"{ref_type}_filtered_reads"]
        deletion_reads = counts[f"{ref_type}_deletion_signal_reads"]
        assessed_bases, gap_bases, gap_total = counts[f"{ref_type}_mpileup_summary"]
        if deletion_reads > filtered_reads:
            raise ValueError(
                f"{args.file_index}/{ref_type}: deletion-bearing reads exceed filtered reads."
            )
        if gap_bases > assessed_bases:
            raise ValueError(
                f"{args.file_index}/{ref_type}: gap-bearing bases exceed assessed bases."
            )
        stats[f"{ref_type}_deletion_signal_reads"] = deletion_reads
        stats[f"{ref_type}_deletion_signal_read_ratio"] = deletion.format_ratio(
            deletion_reads, filtered_reads
        )
        stats[f"{ref_type}_assessed_bases"] = assessed_bases
        stats[f"{ref_type}_gap_signal_bases"] = gap_bases
        stats[f"{ref_type}_gap_signal_total"] = gap_total
        stats[f"{ref_type}_gap_signal_base_ratio"] = deletion.format_ratio(
            gap_bases, assessed_bases
        )

    write_log(args.output, stats)
    print(f"[library-log] Wrote {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, ValueError) as error:
        print(f"[library-log] ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
