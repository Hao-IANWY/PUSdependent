#!/usr/bin/env python3

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Container, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


DEFAULT_KEY_COLUMNS = ("chr", "pos", "strand")


@dataclass(frozen=True)
class SiteRecord:
    key: Tuple[str, ...]
    context: str
    wt_ratio: float


@dataclass
class LineResult:
    name: str
    sites: Dict[Tuple[str, ...], Dict[str, object]]
    all_sites: Dict[Tuple[str, ...], Dict[str, object]]


def parse_file_list(value: str) -> List[Path]:
    files = [Path(item.strip()) for item in value.split(",") if item.strip()]
    if not files:
        raise argparse.ArgumentTypeError("file list cannot be empty")
    return files


def positive_float(value: str) -> float:
    parsed = float(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("value must be non-negative")
    return parsed


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Call dependent sites for one genotype using a WT candidate TSV and two "
            "independent mutant lines, based purely on the modification-ratio drop "
            "effect size (WT_Ratio - Mut_Ratio >= min_delta). No significance test "
            "is performed; a site is dependent when both mutant lines clear the "
            "delta floor and the mock-depth filter."
        )
    )
    parser.add_argument(
        "--tsv",
        required=True,
        type=Path,
        help=(
            "WT candidate TSV with a header. It must contain coordinate columns and "
            "a WT ratio column. Rows in this file are treated as WT candidate sites."
        ),
    )
    parser.add_argument(
        "--genotype1",
        required=True,
        help="Name of the first mutant line, for example trua1a_1-Leaf.",
    )
    parser.add_argument(
        "--genotype1-treat-files",
        required=True,
        type=parse_file_list,
        help="Comma-separated Treat adjusted TSV files for genotype1.",
    )
    parser.add_argument(
        "--genotype1-mock-files",
        required=True,
        type=parse_file_list,
        help="Comma-separated Mock adjusted TSV files for genotype1.",
    )
    parser.add_argument(
        "--genotype2",
        required=True,
        help="Name of the second mutant line, for example trua1a_2-Leaf.",
    )
    parser.add_argument(
        "--genotype2-treat-files",
        required=True,
        type=parse_file_list,
        help="Comma-separated Treat adjusted TSV files for genotype2.",
    )
    parser.add_argument(
        "--genotype2-mock-files",
        required=True,
        type=parse_file_list,
        help="Comma-separated Mock adjusted TSV files for genotype2.",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory for line-level, final-level, and summary files.",
    )
    parser.add_argument(
        "--output-prefix",
        default=None,
        help=(
            "Prefix for final output files. Default: common prefix inferred from "
            "--genotype1 and --genotype2, or genotype1_vs_genotype2."
        ),
    )
    parser.add_argument(
        "--key-cols",
        default=",".join(DEFAULT_KEY_COLUMNS),
        help=(
            "Comma-separated site key columns in adjusted TSV files and output. "
            "Default: chr,pos,strand."
        ),
    )
    parser.add_argument(
        "--tsv-key-cols",
        default=None,
        help=(
            "Comma-separated site key columns in --tsv. Default: use --key-cols. "
            "For transcript-coordinate WT TSVs, use chromosome,position,strand."
        ),
    )
    parser.add_argument(
        "--context-col",
        default="context",
        help="Optional context column name. If absent, context is written as NA.",
    )
    parser.add_argument(
        "--ratio-col",
        default="ratio",
        help="WT ratio column in --tsv. Default: ratio.",
    )
    parser.add_argument(
        "--adjusted-ratio-col",
        default="ratio",
        help="Ratio column in adjusted TSV files. Default: ratio.",
    )
    parser.add_argument(
        "--adjusted-depth-col",
        default="depth",
        help="Depth column in adjusted TSV files. Default: depth.",
    )
    parser.add_argument(
        "--min-mock-depth",
        type=positive_float,
        default=10.0,
        help="Minimum mean Mock depth required per mutant line. Default: 10.",
    )
    parser.add_argument(
        "--min-delta",
        type=positive_float,
        default=0.05,
        help=(
            "Minimum WT_Ratio - Mut_Ratio required per mutant line for a site to be "
            "dependent (effect-size floor, 5%% absolute difference by default). "
            "Default: 0.05."
        ),
    )
    parser.add_argument(
        "--allow-duplicate-sites",
        action="store_true",
        help=(
            "Allow duplicate keys by keeping the first occurrence. Default: fail on "
            "duplicates because duplicate coordinates make joins ambiguous."
        ),
    )
    return parser.parse_args(argv)


def split_columns(value: str, option_name: str) -> Tuple[str, ...]:
    columns = tuple(item.strip() for item in value.split(",") if item.strip())
    if not columns:
        raise ValueError(f"{option_name} must define at least one column")
    return columns


def require_columns(path: Path, fieldnames: Optional[Sequence[str]], columns: Iterable[str]) -> None:
    if fieldnames is None:
        raise ValueError(f"{path}: missing TSV header")
    missing = [column for column in columns if column not in fieldnames]
    if missing:
        raise ValueError(f"{path}: missing required column(s): {', '.join(missing)}")


def column_indices(path: Path, fieldnames: Sequence[str], columns: Iterable[str]) -> Dict[str, int]:
    require_columns(path, fieldnames, columns)
    return {column: fieldnames.index(column) for column in columns}


def read_header(path: Path) -> List[str]:
    with path.open("r", newline="") as handle:
        first_line = handle.readline()
    if not first_line:
        raise ValueError(f"{path}: missing TSV header")
    return first_line.rstrip("\r\n").split("\t")


def to_float(value: str, path: Path, row_number: int, column: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"{path}: row {row_number}: column {column} is not numeric: {value!r}"
        ) from exc


def site_key(row: Mapping[str, str], key_cols: Sequence[str]) -> Tuple[str, ...]:
    return tuple(row[column] for column in key_cols)


def read_wt_sites(
    path: Path,
    key_cols: Sequence[str],
    context_col: str,
    ratio_col: str,
    allow_duplicate_sites: bool,
) -> Dict[Tuple[str, ...], SiteRecord]:
    required_cols = [*key_cols, ratio_col]
    sites: Dict[Tuple[str, ...], SiteRecord] = {}

    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(path, reader.fieldnames, required_cols)
        has_context = context_col in (reader.fieldnames or [])

        for row_number, row in enumerate(reader, start=2):
            key = site_key(row, key_cols)
            if key in sites:
                if allow_duplicate_sites:
                    continue
                raise ValueError(f"{path}: duplicate site key at row {row_number}: {key}")

            context = row[context_col] if has_context else "NA"
            sites[key] = SiteRecord(
                key=key,
                context=context,
                wt_ratio=to_float(row[ratio_col], path, row_number, ratio_col),
            )

    return sites


def read_adjusted_file(
    path: Path,
    key_cols: Sequence[str],
    ratio_col: str,
    depth_col: str,
    allow_duplicate_sites: bool,
    target_keys: Optional[Container[Tuple[str, ...]]] = None,
) -> Dict[Tuple[str, ...], Tuple[float, float]]:
    required_cols = [*key_cols, ratio_col, depth_col]
    values: Dict[Tuple[str, ...], Tuple[float, float]] = {}

    with path.open("r", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        fieldnames = next(reader, None)
        if fieldnames is None:
            raise ValueError(f"{path}: missing TSV header")
        indices = column_indices(path, fieldnames, required_cols)
        key_indices = [indices[column] for column in key_cols]
        ratio_index = indices[ratio_col]
        depth_index = indices[depth_col]

        for row_number, row in enumerate(reader, start=2):
            key = tuple(row[index] for index in key_indices)
            if target_keys is not None and key not in target_keys:
                continue
            if key in values:
                if allow_duplicate_sites:
                    continue
                raise ValueError(f"{path}: duplicate site key at row {row_number}: {key}")

            ratio = to_float(row[ratio_index], path, row_number, ratio_col)
            depth = to_float(row[depth_index], path, row_number, depth_col)
            values[key] = (ratio, depth)
            if target_keys is not None and len(values) == len(target_keys):
                break

    return values


def mean(values: Sequence[float]) -> float:
    if not values:
        raise ValueError("cannot calculate mean of an empty list")
    return sum(values) / len(values)


def call_line_dependent(
    name: str,
    wt_sites: Mapping[Tuple[str, ...], SiteRecord],
    mock_files: Sequence[Path],
    treat_files: Sequence[Path],
    key_cols: Sequence[str],
    adjusted_ratio_col: str,
    adjusted_depth_col: str,
    min_mock_depth: float,
    min_delta: float,
    allow_duplicate_sites: bool,
) -> LineResult:
    target_keys = frozenset(wt_sites)
    mock_tables = [
        read_adjusted_file(
            path,
            key_cols,
            adjusted_ratio_col,
            adjusted_depth_col,
            allow_duplicate_sites,
            target_keys,
        )
        for path in mock_files
    ]
    treat_tables = [
        read_adjusted_file(
            path,
            key_cols,
            adjusted_ratio_col,
            adjusted_depth_col,
            allow_duplicate_sites,
            target_keys,
        )
        for path in treat_files
    ]

    dependent: Dict[Tuple[str, ...], Dict[str, object]] = {}
    all_sites: Dict[Tuple[str, ...], Dict[str, object]] = {}

    for key, wt_site in wt_sites.items():
        if any(key not in table for table in mock_tables):
            continue
        if any(key not in table for table in treat_tables):
            continue

        mock_mean_depth = mean([table[key][1] for table in mock_tables])
        mut_ratios = [table[key][0] for table in treat_tables]
        mut_ratio = mean(mut_ratios)
        delta = wt_site.wt_ratio - mut_ratio

        # 纯效应量判定：不做显著性检验，仅要求修饰比例下降幅度超过 min_delta。
        is_dependent = mock_mean_depth >= min_mock_depth and delta >= min_delta

        site_data = {
            "key": key,
            "context": wt_site.context,
            "WT_Ratio": wt_site.wt_ratio,
            "Mock_Mean_Depth": mock_mean_depth,
            "Mut_Ratio": mut_ratio,
            "Delta_Ratio(WT-Mut)": delta,
            "is_dependent": is_dependent,
        }

        all_sites[key] = site_data
        if is_dependent:
            dependent[key] = site_data

    return LineResult(name=name, sites=dependent, all_sites=all_sites)


def merge_line_results(
    wt_sites: Mapping[Tuple[str, ...], SiteRecord],
    line1: LineResult,
    line2: LineResult,
) -> Dict[Tuple[str, ...], Dict[str, object]]:
    final_sites: Dict[Tuple[str, ...], Dict[str, object]] = {}

    common_keys = line1.all_sites.keys() & line2.all_sites.keys()

    for key in sorted(common_keys):
        wt_site = wt_sites[key]
        is_dep1 = line1.all_sites[key]["is_dependent"]
        is_dep2 = line2.all_sites[key]["is_dependent"]

        mut_ratio = mean(
            [
                float(line1.all_sites[key]["Mut_Ratio"]),
                float(line2.all_sites[key]["Mut_Ratio"]),
            ]
        )
        final_sites[key] = {
            "key": key,
            "context": wt_site.context,
            "WT_Ratio": wt_site.wt_ratio,
            "Mut_Ratio": mut_ratio,
            "Delta_Ratio(WT-Mut)": wt_site.wt_ratio - mut_ratio,
            "is_dependent": (is_dep1 and is_dep2),
        }
    return final_sites


def format_float(value: object) -> str:
    return f"{float(value):.10g}"


def write_line_table(
    path: Path,
    key_cols: Sequence[str],
    line_result: LineResult,
) -> None:
    fieldnames = [
        *key_cols,
        "context",
        "WT_Ratio",
        "Mock_Mean_Depth",
        "Mut_Ratio",
        "Delta_Ratio(WT-Mut)",
        "is_dependent",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(line_result.all_sites):
            row = line_result.all_sites[key]
            writer.writerow(
                {
                    **dict(zip(key_cols, key)),
                    "context": row["context"],
                    "WT_Ratio": format_float(row["WT_Ratio"]),
                    "Mock_Mean_Depth": format_float(row["Mock_Mean_Depth"]),
                    "Mut_Ratio": format_float(row["Mut_Ratio"]),
                    "Delta_Ratio(WT-Mut)": format_float(row["Delta_Ratio(WT-Mut)"]),
                    "is_dependent": str(row["is_dependent"]),
                }
            )


def write_final_table(
    path: Path,
    key_cols: Sequence[str],
    final_sites: Mapping[Tuple[str, ...], Mapping[str, object]],
) -> None:
    fieldnames = [
        *key_cols,
        "context",
        "WT_Ratio",
        "Mut_Ratio",
        "Delta_Ratio(WT-Mut)",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(final_sites):
            row = final_sites[key]
            if not row["is_dependent"]:
                continue
            writer.writerow(
                {
                    **dict(zip(key_cols, key)),
                    "context": row["context"],
                    "WT_Ratio": format_float(row["WT_Ratio"]),
                    "Mut_Ratio": format_float(row["Mut_Ratio"]),
                    "Delta_Ratio(WT-Mut)": format_float(row["Delta_Ratio(WT-Mut)"]),
                }
            )


def infer_output_prefix(genotype1: str, genotype2: str) -> str:
    line_pattern = re.compile(r"^(?P<base>.+?)[_-](?P<line>[12])(?P<suffix>[-_.].*)?$")
    match1 = line_pattern.match(genotype1)
    match2 = line_pattern.match(genotype2)
    if match1 and match2:
        if (
            match1.group("base") == match2.group("base")
            and match1.group("suffix") == match2.group("suffix")
            and match1.group("line") != match2.group("line")
        ):
            return f"{match1.group('base')}{match1.group('suffix') or ''}"

    prefix_chars: List[str] = []
    for left, right in zip(genotype1, genotype2):
        if left != right:
            break
        prefix_chars.append(left)
    prefix = "".join(prefix_chars).rstrip("_-.")
    if prefix:
        return prefix
    return f"{genotype1}_vs_{genotype2}"


def write_summary(
    path: Path,
    args: argparse.Namespace,
    key_cols: Sequence[str],
    tsv_key_cols: Sequence[str],
    wt_count: int,
    line1: LineResult,
    line2: LineResult,
    final_count: int,
    final_all_sites_count: int,
) -> None:
    summary = {
        "tsv": str(args.tsv),
        "key_cols": list(key_cols),
        "tsv_key_cols": list(tsv_key_cols),
        "ratio_col": args.ratio_col,
        "adjusted_ratio_col": args.adjusted_ratio_col,
        "adjusted_depth_col": args.adjusted_depth_col,
        "min_mock_depth": args.min_mock_depth,
        "min_delta": args.min_delta,
        "method": "effect_size",
        "wt_candidate_count": wt_count,
        "genotype1": args.genotype1,
        "genotype1_mock_files": [str(path) for path in args.genotype1_mock_files],
        "genotype1_treat_files": [str(path) for path in args.genotype1_treat_files],
        "genotype1_all_sites_count": len(line1.all_sites),
        "genotype1_dependent_count": len(line1.sites),
        "genotype2": args.genotype2,
        "genotype2_mock_files": [str(path) for path in args.genotype2_mock_files],
        "genotype2_treat_files": [str(path) for path in args.genotype2_treat_files],
        "genotype2_all_sites_count": len(line2.all_sites),
        "genotype2_dependent_count": len(line2.sites),
        "final_all_sites_count": final_all_sites_count,
        "final_dependent_count": final_count,
    }
    with path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run(args: argparse.Namespace) -> int:
    key_cols = split_columns(args.key_cols, "--key-cols")
    tsv_key_cols = split_columns(args.tsv_key_cols or args.key_cols, "--tsv-key-cols")
    if len(tsv_key_cols) != len(key_cols):
        raise ValueError(
            "--tsv-key-cols and --key-cols must define the same number of columns: "
            f"{len(tsv_key_cols)} != {len(key_cols)}"
        )
    args.outdir.mkdir(parents=True, exist_ok=True)

    wt_sites = read_wt_sites(
        args.tsv,
        tsv_key_cols,
        args.context_col,
        args.ratio_col,
        args.allow_duplicate_sites,
    )
    if not wt_sites:
        raise ValueError(f"{args.tsv}: no WT candidate sites found")

    line1 = call_line_dependent(
        args.genotype1,
        wt_sites,
        args.genotype1_mock_files,
        args.genotype1_treat_files,
        key_cols,
        args.adjusted_ratio_col,
        args.adjusted_depth_col,
        args.min_mock_depth,
        args.min_delta,
        args.allow_duplicate_sites,
    )
    line2 = call_line_dependent(
        args.genotype2,
        wt_sites,
        args.genotype2_mock_files,
        args.genotype2_treat_files,
        key_cols,
        args.adjusted_ratio_col,
        args.adjusted_depth_col,
        args.min_mock_depth,
        args.min_delta,
        args.allow_duplicate_sites,
    )
    final_sites = merge_line_results(wt_sites, line1, line2)

    output_prefix = args.output_prefix or infer_output_prefix(args.genotype1, args.genotype2)
    line1_output = args.outdir / f"{args.genotype1}.dependent-sites.tsv"
    line2_output = args.outdir / f"{args.genotype2}.dependent-sites.tsv"
    final_output = args.outdir / f"{output_prefix}.dependent.tsv"
    summary_output = args.outdir / f"{output_prefix}.dependent.summary.json"

    write_line_table(line1_output, key_cols, line1)
    write_line_table(line2_output, key_cols, line2)
    write_final_table(final_output, key_cols, final_sites)
    write_summary(
        summary_output,
        args,
        key_cols,
        tsv_key_cols,
        len(wt_sites),
        line1,
        line2,
        sum(1 for site in final_sites.values() if site["is_dependent"]),
        len(final_sites),
    )

    sys.stderr.write(f"WT candidates: {len(wt_sites)}\n")
    sys.stderr.write(f"{args.genotype1} all sites: {len(line1.all_sites)}\n")
    sys.stderr.write(f"{args.genotype1} dependent sites: {len(line1.sites)}\n")
    sys.stderr.write(f"{args.genotype2} all sites: {len(line2.all_sites)}\n")
    sys.stderr.write(f"{args.genotype2} dependent sites: {len(line2.sites)}\n")
    sys.stderr.write(f"Final all sites: {len(final_sites)}\n")
    sys.stderr.write(f"Final dependent sites: {sum(1 for site in final_sites.values() if site['is_dependent'])}\n")
    sys.stderr.write(f"Final output: {final_output}\n")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
