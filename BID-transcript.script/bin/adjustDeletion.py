#!/home/HL/miniforge3/envs/DeepPSU/bin/python
# -*- coding: utf-8 -*-
#
# Shift deletion signals to the rightmost base of a contiguous T run.

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable, List, Optional, Sequence, Tuple

import pysam


MATCH_OPS = {0, 7, 8}  # M, =, X
REF_CONSUMING_OPS = {0, 2, 3, 7, 8}  # M, D, N, =, X


def positive_int(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "--threads must be a positive integer."
        ) from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("--threads must be a positive integer.")
    return parsed


def normalize_cigartuples(
    cigartuples: Optional[Sequence[Tuple[int, int]]]
) -> List[Tuple[int, int]]:
    if not cigartuples:
        return []

    merged: List[Tuple[int, int]] = []
    for op, length in cigartuples:
        if length <= 0:
            continue
        if merged and merged[-1][0] == op:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    return merged


def cigar_to_string(cigartuples: Sequence[Tuple[int, int]]) -> str:
    op_map = {
        0: "M",
        1: "I",
        2: "D",
        3: "N",
        4: "S",
        5: "H",
        6: "P",
        7: "=",
        8: "X",
    }
    return "".join(f"{length}{op_map[op]}" for op, length in cigartuples)


def build_prefix_ref_positions(
    cigartuples: Sequence[Tuple[int, int]]
) -> List[int]:
    ref_pos = 0
    prefix: List[int] = []
    for op, length in cigartuples:
        prefix.append(ref_pos)
        if op in REF_CONSUMING_OPS:
            ref_pos += length
    return prefix


def get_t_runs(seq: str) -> List[Tuple[int, int]]:
    runs: List[Tuple[int, int]] = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] != "T":
            i += 1
            continue
        j = i + 1
        while j < n and seq[j] == "T":
            j += 1
        runs.append((i, j))
        i = j
    return runs


def interval_has_only_supported_ops(
    cigartuples: Sequence[Tuple[int, int]],
    ref_start: int,
    interval_start: int,
    interval_end: int,
) -> Optional[int]:
    """
    Return the total deletion length inside the interval if the interval is eligible.

    The interval is eligible only when:
    - at least one D overlaps the interval;
    - all overlapping ref-consuming ops are M/=/X or D;
    - no insertion/softclip/hardclip/pad occurs while the reference pointer is inside
      the interval;
    - N is not allowed inside the interval.
    """
    ref_pos = ref_start
    del_total = 0
    has_deletion = False

    for op, length in cigartuples:
        if op in REF_CONSUMING_OPS:
            op_start = ref_pos
            op_end = ref_pos + length
            overlaps = not (op_end <= interval_start or op_start >= interval_end)

            if overlaps:
                if op == 3:
                    return None
                if op == 2:
                    overlap_start = max(op_start, interval_start)
                    overlap_end = min(op_end, interval_end)
                    del_total += overlap_end - overlap_start
                    has_deletion = True
                elif op not in MATCH_OPS:
                    return None

            ref_pos = op_end
        else:
            if interval_start <= ref_pos < interval_end:
                return None

    if not has_deletion:
        return None
    return del_total


def rewrite_interval(
    cigartuples: Sequence[Tuple[int, int]],
    ref_start: int,
    interval_start: int,
    interval_end: int,
    del_total: int,
) -> Optional[List[Tuple[int, int]]]:
    replacement: List[Tuple[int, int]] = []
    match_len = (interval_end - interval_start) - del_total
    if match_len > 0:
        replacement.append((0, match_len))
    replacement.append((2, del_total))

    new_cigar: List[Tuple[int, int]] = []
    ref_pos = ref_start
    inserted = False

    for op, length in cigartuples:
        if op in REF_CONSUMING_OPS:
            op_start = ref_pos
            op_end = ref_pos + length

            if op_end <= interval_start:
                new_cigar.append((op, length))
            elif op_start >= interval_end:
                if not inserted:
                    return None
                new_cigar.append((op, length))
            else:
                if op_start < interval_start:
                    left_len = interval_start - op_start
                    new_cigar.append((op, left_len))
                if not inserted:
                    new_cigar.extend(replacement)
                    inserted = True
                if op_end > interval_end:
                    right_len = op_end - interval_end
                    new_cigar.append((op, right_len))

            ref_pos = op_end
        else:
            if interval_start <= ref_pos < interval_end:
                return None
            new_cigar.append((op, length))

    if not inserted:
        return None
    return normalize_cigartuples(new_cigar)


def shift_read(record: pysam.AlignedSegment, fasta: pysam.FastaFile) -> Tuple[pysam.AlignedSegment, List[dict]]:
    """
    Rewrite a read in place.

    Returns:
        (record, events)
        events is a list of dicts suitable for report TSV output.
    """
    events: List[dict] = []
    if record.is_unmapped or record.reference_id < 0 or not record.cigartuples:
        return record, events

    if record.reference_name not in fasta.references:
        return record, events

    cigar = normalize_cigartuples(record.cigartuples)
    ref_start = record.reference_start
    ref_end = record.reference_end
    if ref_start is None or ref_end is None or ref_end <= ref_start:
        return record, events

    ref_seq = fasta.fetch(record.reference_name, ref_start, ref_end).upper()
    t_runs = get_t_runs(ref_seq)
    if not t_runs:
        return record, events

    for run_rel_start, run_rel_end in t_runs:
        interval_start = ref_start + run_rel_start
        interval_end = ref_start + run_rel_end

        del_total = interval_has_only_supported_ops(
            cigar, ref_start, interval_start, interval_end
        )
        if del_total is None:
            continue

        new_cigar = rewrite_interval(
            cigar, ref_start, interval_start, interval_end, del_total
        )
        if new_cigar is None:
            continue
        if new_cigar == cigar:
            continue

        old_cigar = cigar_to_string(cigar)
        cigar = new_cigar
        new_cigar_str = cigar_to_string(cigar)
        events.append(
            {
                "read_name": record.query_name,
                "chrom": record.reference_name,
                "interval_start": interval_start + 1,
                "interval_end": interval_end,
                "old_cigar": old_cigar,
                "new_cigar": new_cigar_str,
                "deleted_bases": del_total,
            }
        )

    record.cigartuples = cigar
    return record, events


def open_report(report_path: Optional[str]):
    if report_path is None:
        return None
    os.makedirs(os.path.dirname(os.path.abspath(report_path)) or ".", exist_ok=True)
    report = open(report_path, "w", encoding="utf-8")
    report.write(
        "read_name\tchrom\tinterval_start_1based\tinterval_end_1based\told_cigar\tnew_cigar\tdeleted_bases\n"
    )
    return report


def write_report(report, events: Iterable[dict]) -> None:
    if report is None:
        return
    for event in events:
        report.write(
            "{read_name}\t{chrom}\t{interval_start}\t{interval_end}\t{old_cigar}\t{new_cigar}\t{deleted_bases}\n".format(
                **event
            )
        )


def run(
    input_bam: str,
    reference: str,
    output_bam: str,
    report_tsv: Optional[str],
    threads: int = 1,
) -> None:
    samtools = shutil.which("samtools")
    if samtools is None:
        raise SystemExit("samtools not found in PATH.")
    if threads < 1:
        raise SystemExit("--threads must be a positive integer.")

    os.makedirs(os.path.dirname(os.path.abspath(output_bam)) or ".", exist_ok=True)
    fasta = pysam.FastaFile(reference)
    temp_dir = os.path.dirname(os.path.abspath(output_bam)) or "."
    temp_bam = tempfile.NamedTemporaryFile(
        prefix="adjustDeletion_", suffix=".bam", dir=temp_dir, delete=False
    )
    temp_bam_path = temp_bam.name
    temp_bam.close()

    report = open_report(report_tsv)
    total_reads = 0
    shifted_reads = 0
    shifted_events = 0

    try:
        with pysam.AlignmentFile(
            input_bam, "rb", threads=threads
        ) as fin, pysam.AlignmentFile(
            temp_bam_path, "wb", template=fin, threads=threads
        ) as fout:
            for record in fin:
                total_reads += 1
                new_record, events = shift_read(record, fasta)
                if events:
                    shifted_reads += 1
                    shifted_events += len(events)
                    write_report(report, events)
                fout.write(new_record)

        if report is not None:
            report.close()

        with open(output_bam, "wb") as out_fh:
            subprocess.run(
                [samtools, "calmd", "-@", str(threads), "-b", temp_bam_path, reference],
                stdout=out_fh,
                stderr=subprocess.DEVNULL,
                check=True,
            )

        subprocess.run([samtools, "index", "-@", str(threads), output_bam], check=True)
        sys.stderr.write(
            f"Processed {total_reads} reads, shifted {shifted_events} deletion event(s) in {shifted_reads} read(s).\n"
        )
    finally:
        try:
            os.remove(temp_bam_path)
        except FileNotFoundError:
            pass
        if report is not None and not report.closed:
            report.close()


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Shift deletion signals to the rightmost base of a contiguous T run in BAM."
        )
    )
    parser.add_argument(
        "-i",
        "--input-bam",
        required=True,
        help="Input BAM file.",
    )
    parser.add_argument(
        "-r",
        "--reference",
        required=True,
        help="Reference FASTA file.",
    )
    parser.add_argument(
        "-o",
        "--output-bam",
        required=True,
        help="Output BAM file.",
    )
    parser.add_argument(
        "--report-tsv",
        default=None,
        help="Optional TSV report of shifted deletion events.",
    )
    parser.add_argument(
        "--threads",
        type=positive_int,
        default=1,
        help="Number of threads to use for BAM I/O, calmd, and index.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    run(
        args.input_bam,
        args.reference,
        args.output_bam,
        args.report_tsv,
        args.threads,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
