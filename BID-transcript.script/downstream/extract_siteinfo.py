#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extracts specific site attributes (depth, gap, etc.) from multiple adjusted.tsv files.

Background
----------
`STAT/adjustGap/<ref_type>/<sample>.adjusted.tsv` files are transcriptome-wide 
base-resolution tables (transcript ID, 1-based coordinate, >40M rows per file).
This tool performs precise (chr, pos) matching to aggregate per-sample strands, 
contexts, ratios, depths, and gaps into a single summary table. Missing data is 
reported as NA. Hits per file are reported to stderr.

Input is a YAML configuration:
    sites:
      - chr: AT1G06190.5
        pos: 48
      - AT1G06190.5:66
    files:
      - Temperature/STAT/adjustGap/mRNA/Mock-21-trua1b_1.rep1.adjusted.tsv
      - Temperature/STAT/adjustGap/mRNA/Mock-27-*.rep?.adjusted.tsv
    format: wide                # Optional: wide (default) / long
    output: siteinfo.tsv        # Optional; defaults to stdout

Output Formats
--------------
wide — One site per row: `chr pos strand context` followed by 
       `<sample>.ratio` / `<sample>.depth` / `<sample>.gap` columns for each sample.
long — One site × sample combination per row: `chr pos sample strand context ratio depth gap`.

Relative paths are resolved relative to the YAML file, then the repository root.

Usage
-----
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

# adjusted.tsv columns: chr, pos, strand, context, ratio, depth, gap
FIELDS = ("strand", "context", "ratio", "depth", "gap")
# Per-sample attributes in wide format
PER_SAMPLE = ("ratio", "depth", "gap")
MISSING = ("NA",) * len(FIELDS)


def parse_site_entry(entry: object) -> Tuple[str, int]:
    """Parse site entry from YAML into (chr, pos)."""
    if isinstance(entry, str):
        # Format "chr:pos"
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
    """Resolve file entries to absolute paths."""
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
    """Scan file using grep -F and return site data."""
    argv = ["grep", "-F"]
    for chr_id, pos in queries:
        argv += ["-e", f"{chr_id}\t{pos}\t"]
    argv.append(str(path))
    proc = subprocess.run(argv, capture_output=True, text=True, env={**os.environ, "LC_ALL": "C"})
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"grep failed on {path} (exit {proc.returncode}): {proc.stderr.strip()}")

    hits: Dict[Tuple[str, int], Tuple[str, ...]] = {}
    query_keys = set(queries)
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        key = (fields[0], int(fields[1]))
        if key not in query_keys:
            continue
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
    parser.add_argument("config", type=Path, help="YAML configuration file.")
    parser.add_argument("-o", "--output", type=Path, help="Output TSV file.")
    args = parser.parse_args(argv)

    with open(args.config) as handle:
        config = safe_load(handle) or {}
    for key in ("sites", "files"):
        if not config.get(key):
            raise ValueError(f"missing or empty '{key}' in {args.config}")

    queries = [parse_site_entry(entry) for entry in config["sites"]]
    files = resolve_input_files(config["files"], (args.config.resolve().parent, REPO_ROOT))

    name_counts = Counter(sample_name(path) for path in files)
    per_file = []
    for path in files:
        sample = sample_name(path)
        if name_counts[sample] > 1:
            sample = f"{path.parent.name}.{sample}"
        hits = extract_sites(path, queries)
        per_file.append((sample, hits))
        print(f"{path.name}: {sum(1 for q in queries if q in hits)}/{len(queries)} sites found",
              file=sys.stderr)

    if config.get("format", "wide") == "long":
        header = ["chr", "pos", "sample", *FIELDS]
        rows = [[chr_id, pos, sample, *hits.get((chr_id, pos), MISSING)]
                for chr_id, pos in queries for sample, hits in per_file]
    else:
        header = ["chr", "pos", "strand", "context"]
        header += [f"{sample}.{key}" for sample, _ in per_file for key in PER_SAMPLE]
        rows = []
        for chr_id, pos in queries:
            site = (chr_id, pos)
            ref = next((hits[site][:2] for _, hits in per_file if site in hits), ("NA", "NA"))
            row = [chr_id, pos, *ref]
            for _, hits in per_file:
                fields = hits.get(site)
                row += list(fields[2:5]) if fields else ["NA"] * len(PER_SAMPLE)
            rows.append(row)

    if args.output:
        output = args.output
    else:
        output = config.get("output")
        if output and not Path(output).is_absolute():
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
