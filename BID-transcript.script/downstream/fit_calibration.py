#!/usr/bin/env python3

import argparse
import csv
import gzip
import math
import sys
from pathlib import Path
from typing import Dict, Mapping, Optional, Sequence, TextIO, Tuple


CurveParams = Tuple[float, float, float]
# Calibration_curve.tsv 与本脚本放在同一目录，按脚本位置解析以避免硬编码绝对路径。
DEFAULT_CURVE = Path(__file__).resolve().parent / "Calibration_curve.tsv"


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Convert gap ratio values to calibrated fraction values using "
            "motif-specific c/s/b calibration parameters."
        )
    )
    parser.add_argument(
        "--curve",
        default=DEFAULT_CURVE,
        type=Path,
        help=(
            "Calibration TSV with Motif, c, s, and b columns. "
            f"Default: {DEFAULT_CURVE}"
        ),
    )
    parser.add_argument(
        "--tsv",
        required=True,
        type=Path,
        help="Input TSV with a header containing context and ratio columns.",
    )
    parser.add_argument(
        "--output",
        required=True,
        type=Path,
        help="Output TSV path. Use '-' to write to stdout.",
    )
    parser.add_argument(
        "--fraction-col",
        default="fraction",
        help="Output fraction column name. Default: fraction.",
    )
    parser.add_argument(
        "--digits",
        type=int,
        default=3,
        help="Number of decimal digits for finite fraction values. Default: 3.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail when a row context is absent from the calibration curve.",
    )
    return parser.parse_args(argv)


def normalize_motif(value: str) -> str:
    return value.strip().upper().replace("U", "T").replace("-", "T")


def open_input(path: Path) -> TextIO:
    if str(path) == "-":
        return sys.stdin
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def open_output(path: Path) -> TextIO:
    if str(path) == "-":
        return sys.stdout
    if path.suffix == ".gz":
        return gzip.open(path, "wt", newline="")
    return path.open("w", newline="")


def require_columns(
    path: Path, fieldnames: Optional[Sequence[str]], columns: Sequence[str]
) -> None:
    if fieldnames is None:
        raise ValueError(f"{path}: missing TSV header")
    missing = [column for column in columns if column not in fieldnames]
    if missing:
        raise ValueError(f"{path}: missing required column(s): {', '.join(missing)}")


def parse_float(value: str, path: Path, row_number: int, column: str) -> float:
    if value is None or value.strip() == "":
        return math.nan
    try:
        return float(value)
    except ValueError as exc:
        raise ValueError(
            f"{path}: row {row_number}: column {column} is not numeric: {value!r}"
        ) from exc


def read_curve(path: Path) -> Dict[str, CurveParams]:
    params: Dict[str, CurveParams] = {}
    with open_input(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames, ("Motif", "c", "s", "b"))

        for row_number, row in enumerate(reader, start=2):
            motif = normalize_motif(row["Motif"])
            if not motif:
                raise ValueError(f"{path}: row {row_number}: empty Motif value")
            if motif in params:
                raise ValueError(f"{path}: row {row_number}: duplicate Motif: {motif}")
            params[motif] = (
                parse_float(row["c"], path, row_number, "c"),
                parse_float(row["s"], path, row_number, "s"),
                parse_float(row["b"], path, row_number, "b"),
            )
    return params


def fit_calibration(ratio: float, motif: str, params: Mapping[str, CurveParams]) -> float:
    motif = normalize_motif(motif)
    if math.isnan(ratio) or motif not in params:
        return math.nan

    c, s, b = params[motif]
    if ratio <= b:
        return 0.0

    max_ratio = (c * (1 - s) + (1 - c) * b) / (1 - c * s)
    if ratio >= max_ratio:
        return 1.0

    return (b - ratio) / (c * (b + s - s * ratio - 1))


def format_fraction(value: float, digits: int) -> str:
    if math.isnan(value):
        return "NaN"
    return f"{value:.{digits}f}"


def add_fraction(
    tsv_path: Path,
    output_path: Path,
    params: Mapping[str, CurveParams],
    fraction_col: str,
    digits: int,
    strict: bool,
) -> None:
    with open_input(tsv_path) as input_handle, open_output(output_path) as output_handle:
        reader = csv.DictReader(input_handle, delimiter="\t")
        require_columns(tsv_path, reader.fieldnames, ("context", "ratio"))
        fieldnames = list(reader.fieldnames or [])
        if fraction_col not in fieldnames:
            fieldnames.append(fraction_col)

        writer = csv.DictWriter(
            output_handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        for row_number, row in enumerate(reader, start=2):
            context = normalize_motif(row["context"])
            if strict and context not in params:
                raise ValueError(
                    f"{tsv_path}: row {row_number}: context not found in curve: "
                    f"{row['context']!r}"
                )
            ratio = parse_float(row["ratio"], tsv_path, row_number, "ratio")
            row[fraction_col] = format_fraction(
                fit_calibration(ratio, context, params), digits
            )
            writer.writerow(row)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.digits < 0:
        raise ValueError("--digits must be non-negative")

    params = read_curve(args.curve)
    add_fraction(
        tsv_path=args.tsv,
        output_path=args.output,
        params=params,
        fraction_col=args.fraction_col,
        digits=args.digits,
        strict=args.strict,
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"fit_calibration.py: error: {exc}", file=sys.stderr)
        raise SystemExit(1)
