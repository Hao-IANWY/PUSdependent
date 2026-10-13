#!/usr/bin/env python3
"""Concatenate one sample's ordered adjustGap chunk outputs."""

import argparse
import os
import shutil
from pathlib import Path
from typing import Sequence

EXPECTED_HEADER = "chr\tpos\tstrand\tcontext\tratio\tdepth\tgap\n"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--sample",
        required=True,
        help="Sample/file-index name used inside each adjusted chunk directory",
    )
    parser.add_argument(
        "chunk_dirs",
        nargs="+",
        help="Adjusted chunk directories in reference/chunk order",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    sample = Path(args.sample).name
    inputs: Sequence[Path] = [
        Path(chunk_dir) / f"{sample}.adjusted.tsv" for chunk_dir in args.chunk_dirs
    ]
    if not inputs:
        raise ValueError("No adjusted chunk directories supplied")

    missing = [str(path) for path in inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing adjusted chunk input(s): " + ", ".join(missing))

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temp_output = output.with_name(f".{output.name}.tmp.{os.getpid()}")

    try:
        with temp_output.open("w", encoding="utf-8", newline="") as output_handle:
            output_handle.write(EXPECTED_HEADER)
            for path in inputs:
                with path.open("r", encoding="utf-8", newline="") as input_handle:
                    header = input_handle.readline()
                    if header != EXPECTED_HEADER:
                        raise ValueError(f"{path}: unexpected adjustGap header: {header!r}")
                    shutil.copyfileobj(input_handle, output_handle)
        os.replace(temp_output, output)
    except BaseException:
        temp_output.unlink(missing_ok=True)
        raise


if __name__ == "__main__":
    main()
