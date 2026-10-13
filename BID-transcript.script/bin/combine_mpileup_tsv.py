#!/home/HL/miniforge3/envs/DeepPSU/bin/python

import argparse
from pathlib import Path

import polars as pl

KEY_COLUMNS = ["chr", "pos", "ref_base", "strand"]
VALUE_COLUMN = "value"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--suffix", required=True)
    parser.add_argument("inputs", nargs="+")
    return parser.parse_args()


def read_sample_table(file_path: Path, sample_name: str) -> pl.DataFrame:
    return pl.read_csv(
        file_path,
        separator="\t",
        has_header=False,
        new_columns=[*KEY_COLUMNS, VALUE_COLUMN],
    ).rename({VALUE_COLUMN: sample_name})


def main():
    args = parse_args()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    input_files = sorted(Path(file_path) for file_path in args.inputs)
    if not input_files:
        raise FileNotFoundError("未找到可合并的 mpileup tsv 文件")

    first_file = input_files[0]
    if not first_file.name.endswith(args.suffix):
        raise ValueError(f"文件名与 suffix 不匹配: {first_file}")

    first_sample = first_file.name.removesuffix(args.suffix)
    combined_df = read_sample_table(first_file, first_sample)

    for file_path in input_files[1:]:
        if not file_path.name.endswith(args.suffix):
            raise ValueError(f"文件名与 suffix 不匹配: {file_path}")
        sample_name = file_path.name.removesuffix(args.suffix)
        combined_df = combined_df.join(
            read_sample_table(file_path, sample_name),
            on=KEY_COLUMNS,
            how="inner",
        )

    combined_df.write_csv(output_path, separator="\t")


if __name__ == "__main__":
    main()
