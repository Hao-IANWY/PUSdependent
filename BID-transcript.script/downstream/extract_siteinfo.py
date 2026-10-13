#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从一批 adjusted.tsv 中提取指定位点的 depth/gap 等属性。

背景
----
`STAT/adjustGap/<ref_type>/<sample>.adjusted.tsv` 是全转录组逐碱基表
（chr=转录本 ID，pos=转录本内 1-based 坐标，单文件四千万行以上），
人工核对某几个位点时不可能整表打开。本工具按 (chr, pos) 精确匹配，
把每个位点在各输入样品上的 strand/context/ratio/depth/gap 汇总成一张表，
查不到的组合记 NA，stderr 按文件报告命中数。

输入是 yaml，`sites` 支持两种等价写法，`files` 支持 glob 通配符：

    sites:
      - chr: AT1G06190.5          # 详细写法
        pos: 48
      - AT1G06190.5:66            # 简写 "chr:pos"
    files:
      - Temperature/STAT/adjustGap/mRNA/Mock-21-trua1b_1.rep1.adjusted.tsv
      - Temperature/STAT/adjustGap/mRNA/Mock-27-*.rep?.adjusted.tsv
    format: wide                  # 可选：wide（默认）/ long
    output: siteinfo.tsv          # 可选；缺省打印到 stdout

输出两种版式，`format` 控制：
    wide — 一个位点一行：`chr pos strand context` 之后每个样品展开
           `<sample>.ratio` / `<sample>.depth` / `<sample>.gap` 三列；
    long — 一个 site × sample 组合一行：`chr pos sample strand context ratio depth gap`。

相对路径先按 yaml 所在目录解析、再按仓库根目录解析；
文件可混跨 ref_type（mRNA/rRNA/ncRNA 的 chr 集不同，不在该文件里的位点记 NA）。
样品名取自文件名去掉 `.adjusted.tsv` 的部分；输入里出现同名样品时（跨 ref_type 混合）
用上级目录名做前缀消歧（`mRNA.Mock-21-WT.rep1` / `rRNA.Mock-21-WT.rep1`）。

用法
----
    python3 extract_siteinfo.py sites.yaml [-o out.tsv]
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from yaml import safe_load

REPO_ROOT = Path(__file__).resolve().parents[2]
ADJUSTED_SUFFIX = ".adjusted.tsv"
GLOB_CHARS = set("*?[")

# adjusted.tsv 列：chr pos strand context ratio depth gap
FIELDS = ("strand", "context", "ratio", "depth", "gap")
# wide 版式里每个样品展开的属性列（strand/context 只依赖参考序列，提到位点级）
PER_SAMPLE = ("ratio", "depth", "gap")
MISSING = ("NA",) * len(FIELDS)


def parse_site_entry(entry: object) -> Tuple[str, int]:
    """把 yaml 里一条位点解析成 (chr, pos)。"""
    if isinstance(entry, str):
        # 简写 "chr:pos"；转录本 ID 自身不含冒号
        chr_id, sep, pos = entry.partition(":")
        if not sep:
            raise ValueError(f"invalid site '{entry}', expected 'chr:pos' or {{chr, pos}}")
    elif isinstance(entry, dict):
        chr_id = str(entry["chr"])
        pos = str(entry["pos"])
    else:
        raise ValueError(f"invalid site entry: {entry!r}")
    return chr_id, int(pos)


def resolve_input_files(entries: List[str], bases: Tuple[Path, ...]) -> List[Path]:
    """把 yaml 里的文件条目解析成绝对路径。

    含通配符的条目对每个 base 做 glob 展开；普通条目按 base 顺序找第一个存在的。
    两种情况都要么命中、要么报错，不留静默空集。
    """
    resolved: List[Path] = []
    for entry in entries:
        if not isinstance(entry, str) or not entry.strip():
            raise ValueError(f"invalid file entry: {entry!r}")
        if GLOB_CHARS & set(entry):
            hits = list(dict.fromkeys(p.resolve() for base in bases for p in base.glob(entry)))
            if not hits:
                raise FileNotFoundError(f"glob '{entry}' matched nothing under {[str(b) for b in bases]}")
            resolved.extend(sorted(hits))
        else:
            for base in bases:
                candidate = (base / entry).resolve()
                if candidate.exists():
                    resolved.append(candidate)
                    break
            else:
                raise FileNotFoundError(f"file not found: {entry} (searched {[str(b) for b in bases]})")
    return resolved


def extract_sites(path: Path, queries: List[Tuple[str, int]]) -> Dict[Tuple[str, int], Tuple[str, ...]]:
    """一次 grep -F 扫描，返回 {(chr, pos): (strand, context, ratio, depth, gap)}。

    四千万行的表逐行读 Python 要一分钟起步，而 grep 固定串扫描不到一秒；
    模式 `chr\tpos\t` 带两道制表符，天然不会误配 pos 前缀（如 123 命中 1234）。
    """
    argv = ["grep", "-F"]
    for chr_id, pos in queries:
        argv += ["-e", f"{chr_id}\t{pos}\t"]
    argv.append(str(path))
    proc = subprocess.run(argv, capture_output=True, text=True, env={**os.environ, "LC_ALL": "C"})
    if proc.returncode not in (0, 1):  # 1 = 一条都没匹配到，属于正常输入
        raise RuntimeError(f"grep failed on {path} (exit {proc.returncode}): {proc.stderr.strip()}")

    hits: Dict[Tuple[str, int], Tuple[str, ...]] = {}
    query_keys = set(queries)
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        key = (fields[0], int(fields[1]))
        if key not in query_keys:
            continue  # 防御：模式是子串匹配，逐条回验命中行确属查询集
        hits[key] = tuple(fields[2:7])
    return hits


def sample_name(path: Path) -> str:
    if path.name.endswith(ADJUSTED_SUFFIX):
        return path.name[: -len(ADJUSTED_SUFFIX)]
    return path.name


def main(argv: Optional[Sequence[str]] = None) -> None:
    parser = argparse.ArgumentParser(
        description="Extract depth/gap attributes of given sites from one or more adjusted.tsv files."
    )
    parser.add_argument("config", type=Path, help="yaml config with 'sites' and 'files' (see module docstring)")
    parser.add_argument("-o", "--output", type=Path, help="output TSV; overrides 'output' in the yaml")
    args = parser.parse_args(argv)

    with open(args.config) as handle:
        config = safe_load(handle) or {}
    for key in ("sites", "files"):
        if not config.get(key):
            raise ValueError(f"missing or empty '{key}' in {args.config}")

    queries = [parse_site_entry(entry) for entry in config["sites"]]
    files = resolve_input_files(config["files"], (args.config.resolve().parent, REPO_ROOT))

    name_counts = Counter(sample_name(path) for path in files)
    per_file = []  # [(sample, {site: fields}), ...]，按输入文件顺序
    for path in files:
        sample = sample_name(path)
        if name_counts[sample] > 1:  # mRNA/ 与 rRNA/ 下同一样品名各有一份表，用上级目录名消歧
            sample = f"{path.parent.name}.{sample}"
        hits = extract_sites(path, queries)
        per_file.append((sample, hits))
        print(f"{path.name}: {sum(1 for q in queries if q in hits)}/{len(queries)} sites found",
              file=sys.stderr)

    if config.get("format", "wide") == "long":
        header = ["chr", "pos", "sample", *FIELDS]
        rows = [[chr_id, pos, sample, *hits.get((chr_id, pos), MISSING)]
                for chr_id, pos in queries for sample, hits in per_file]
    else:  # wide：一个位点一行，每个样品展开 ratio/depth/gap 三列
        header = ["chr", "pos", "strand", "context"]
        header += [f"{sample}.{key}" for sample, _ in per_file for key in PER_SAMPLE]
        rows = []
        for chr_id, pos in queries:
            site = (chr_id, pos)
            # strand/context 只依赖参考序列，取任一命中样品的值即可
            ref = next((hits[site][:2] for _, hits in per_file if site in hits), ("NA", "NA"))
            row = [chr_id, pos, *ref]
            for _, hits in per_file:
                fields = hits.get(site)
                row += list(fields[2:5]) if fields else ["NA"] * len(PER_SAMPLE)
            rows.append(row)

    if args.output:
        output = args.output  # CLI 给的路径按调用方 CWD 解释
    else:
        output = config.get("output")
        if output and not Path(output).is_absolute():  # yaml 里的 output 按 yaml 所在目录解释
            output = args.config.resolve().parent / output
    handle = open(output, "w", newline="") if output else sys.stdout
    try:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    finally:
        if handle is not sys.stdout:
            handle.close()
    print(f"{len(rows)} rows -> {output if output else 'stdout'}", file=sys.stderr)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, FileNotFoundError, RuntimeError) as exc:
        sys.exit(f"extract_siteinfo.py: error: {exc}")
