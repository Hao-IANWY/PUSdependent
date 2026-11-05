#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
What-if check for the WT landscape positive-site calls:
add a per-replicate Mock gap-ratio ceiling (< 5% in EVERY Mock replicate)
on top of the current per-replicate criteria
(Treat gap >= 5, Treat ratio >= 5%, Treat ratio >= 2 x Mock ratio).

Reads the WT Fisher tables (which carry per-replicate Mock ratios and the
Pass flag) for the three datasets used by Figures/wildtype_landscape:
  polyARNA/STAT/fisher/mRNA/WT-{Root,leaf,Flower}
  polyARNA/STAT/fisher/rRNA/WT-{Root,leaf,Flower}
  tRNA/STAT/fisher/ncRNA/WT_{Root,Leaf,Flower}   (nuclear tRNA only: ^tRNA-)

Outputs a per-dataset summary TSV, a failing-sites TSV, and prints a
per-RNA-type union summary.
"""
import argparse
import os
import sys

import pandas as pd

# The 9 WT datasets behind the wildtype_landscape figures
DATASETS = [
    ("mRNA",  "Root",   "polyARNA/STAT/fisher/mRNA/WT-Root_sites.Fisher.tsv",    None),
    ("mRNA",  "Leaf",   "polyARNA/STAT/fisher/mRNA/WT-leaf_sites.Fisher.tsv",    None),
    ("mRNA",  "Flower", "polyARNA/STAT/fisher/mRNA/WT-Flower_sites.Fisher.tsv",  None),
    ("rRNA",  "Root",   "polyARNA/STAT/fisher/rRNA/WT-Root_sites.Fisher.tsv",    None),
    ("rRNA",  "Leaf",   "polyARNA/STAT/fisher/rRNA/WT-leaf_sites.Fisher.tsv",    None),
    ("rRNA",  "Flower", "polyARNA/STAT/fisher/rRNA/WT-Flower_sites.Fisher.tsv",  None),
    ("tRNA",  "Root",   "tRNA/STAT/fisher/ncRNA/WT_Root_sites.Fisher.tsv",       r"^tRNA-"),
    ("tRNA",  "Leaf",   "tRNA/STAT/fisher/ncRNA/WT_Leaf_sites.Fisher.tsv",       r"^tRNA-"),
    ("tRNA",  "Flower", "tRNA/STAT/fisher/ncRNA/WT_Flower_sites.Fisher.tsv",     r"^tRNA-"),
]

MOCK_CEILING = 0.05  # proposed rule: Mock gap ratio must be below 5% in every replicate


def analyze_dataset(repo_root, rna, tissue, rel_path, chr_regex):
    path = os.path.join(repo_root, rel_path)
    df = pd.read_csv(path, sep="\t", dtype={"chr": str})
    if chr_regex is not None:
        df = df[df["chr"].str.match(chr_regex)].copy()

    mock_ratio_cols = [c for c in df.columns if c.startswith("Mock-") and c.endswith("_ratio")]
    mock_depth_cols = [c for c in df.columns if c.startswith("Mock-") and c.endswith("_depth")]
    mock_ratio_cols.sort(key=lambda c: c)
    mock_depth_cols.sort(key=lambda c: c)

    pos = df[df["Pass"] == 1].copy()
    n_pos = len(pos)

    mr = pos[mock_ratio_cols]
    md = pos[mock_depth_cols]

    max_mock_ratio = mr.max(axis=1)
    n_bad_reps = (mr >= MOCK_CEILING).sum(axis=1)
    fails = n_bad_reps >= 1

    # how bad are the failing Mock replicates: covered (depth >= 10) or thin?
    bad_cell = mr >= MOCK_CEILING
    bad_depth = md.where(bad_cell)
    fail_with_cov = (bad_depth.max(axis=1) >= 10)
    fail_thin = fails & ~fail_with_cov.fillna(False)

    keep = pos[~fails]
    drop = pos[fails]

    summary = {
        "RNA": rna,
        "Tissue": tissue,
        "n_positive": n_pos,
        "n_fail_mock5pct": int(fails.sum()),
        "n_keep": len(keep),
        "pct_dropped": round(100.0 * fails.sum() / n_pos, 2) if n_pos else 0.0,
        "n_fail_mockdepth_ge10": int(fail_with_cov.sum()),
        "n_fail_mockdepth_lt10": int(fail_thin.sum()),
        "n_fail_1rep": int((n_bad_reps == 1).sum()),
        "n_fail_2rep": int((n_bad_reps == 2).sum()),
        "n_fail_3rep": int((n_bad_reps == 3).sum()),
        "n_fail_also_avg_ge5pct": int((pos["Mock_avg_ratio"] >= MOCK_CEILING)[fails].sum()),
        "median_max_mock_ratio_of_fail": round(float(max_mock_ratio[fails].median()), 4) if fails.sum() else None,
    }
    return summary, pos, drop


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-root", required=True,
                    help="PUSdependent repo root (paths of DATASETS are relative to it)")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    repo_root = os.path.abspath(args.repo_root)
    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    summaries, failing_frames = [], []
    for rna, tissue, rel_path, regex in DATASETS:
        summary, pos, drop = analyze_dataset(repo_root, rna, tissue, rel_path, regex)
        summaries.append(summary)
        if len(drop):
            d = drop[["chr", "pos", "strand"] + [c for c in drop.columns if c.startswith("Mock-") and c.endswith("_ratio")]].copy()
            d.insert(0, "Tissue", tissue)
            d.insert(0, "RNA", rna)
            failing_frames.append(d)
        print(f"[{rna}/{tissue}] positive={summary['n_positive']:>6}  "
              f"fail(any Mock rep >= 5%)={summary['n_fail_mock5pct']:>5} "
              f"({summary['pct_dropped']}%)  keep={summary['n_keep']:>6}", flush=True)

    summary_df = pd.DataFrame(summaries)
    summary_path = os.path.join(outdir, "mock_gap5pct_summary.tsv")
    summary_df.to_csv(summary_path, sep="\t", index=False)

    fail_df = pd.concat(failing_frames, ignore_index=True) if failing_frames else pd.DataFrame()
    fail_path = os.path.join(outdir, "mock_gap5pct_dropped_sites.tsv")
    fail_df.to_csv(fail_path, sep="\t", index=False)

    # union across tissues per RNA type (site-level identity: chr, pos, strand)
    print("\n=== Union across tissues (per RNA type) ===")
    union_rows = []
    for rna in ["mRNA", "rRNA", "tRNA"]:
        sets = {t: set() for t in ["cur", "new"]}
        for rna2, tissue, rel_path, regex in DATASETS:
            if rna2 != rna:
                continue
            _, pos, drop = analyze_dataset(repo_root, rna2, tissue, rel_path, regex)
            key_cur = set(map(tuple, pos[["chr", "pos", "strand"]].values))
            key_drop = set(map(tuple, drop[["chr", "pos", "strand"]].values)) if len(drop) else set()
            sets["cur"] |= key_cur
            sets["new"] |= (key_cur - key_drop)
        union_rows.append({
            "RNA": rna,
            "union_current": len(sets["cur"]),
            "union_dropped": len(sets["cur"] - sets["new"]),
            "union_keep": len(sets["new"]),
            "pct_dropped": round(100.0 * (len(sets["cur"]) - len(sets["new"])) / len(sets["cur"]), 2) if sets["cur"] else 0.0,
        })
        print(f"[{rna}] union current={len(sets['cur'])}  dropped={len(sets['cur']) - len(sets['new'])}  "
              f"keep={len(sets['new'])} ({union_rows[-1]['pct_dropped']}% dropped)")

    total_cur = sum(r["union_current"] for r in union_rows)
    total_drop = sum(r["union_dropped"] for r in union_rows)
    print(f"[ALL] union current={total_cur}  dropped={total_drop}  keep={total_cur - total_drop} "
          f"({round(100.0 * total_drop / total_cur, 2)}% dropped)")

    union_df = pd.DataFrame(union_rows)
    union_path = os.path.join(outdir, "mock_gap5pct_union.tsv")
    union_df.to_csv(union_path, sep="\t", index=False)

    print(f"\nWritten:\n  {summary_path}\n  {fail_path}\n  {union_path}")


if __name__ == "__main__":
    sys.exit(main())
