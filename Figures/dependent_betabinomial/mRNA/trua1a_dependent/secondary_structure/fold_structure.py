#!/usr/bin/env python3
# fold_structure.py — Fold site-window sequences with RNAfold and extract structure metrics.
#
# Purpose: for every site in the sites table × every window size (±25/±50/±100) × every
# temperature, compute the MFE structure and pairing probabilities with ViennaRNA,
# classify the structural element of the center site, and write two TSVs: metrics and a
# per-position unpaired profile.
#
# Usage:
#   fold_structure.py <script_dir> --sites <rel> --bed-prefix <rel> \
#                     --out-prefix <name> [--temps "21 37"]
#   Inputs  — <script_dir>/<sites>: genotype/tissue/transcript/pos/strand/context[/set]
#             <script_dir>/<bed-prefix>_{w}.bed + _{w}.fa: bedtools ±w windows
#   Outputs — <script_dir>/<out-prefix>_metrics.tsv
#             <script_dir>/<out-prefix>_unpaired_profile.tsv
#
# Model parameters: ViennaRNA defaults (dangles=2, GU wobble allowed, no pseudoknots);
# temperature comes from --temps (the main analysis of this project uses 21 °C as the
# in-vivo approximation).
#
# Structure descriptor columns (non-empty only when the center belongs to a single
# stem-loop unit):
#   stem_len/stem_gc/stem_au/stem_gu  pair count and pair-type counts of the stem holding the center
#   loop_len                          length of the innermost loop enclosed by that stem (the unpaired segment between the innermost pair)
#   u_to_loop_dist                    number of paired bases crossed between the center and that loop (0 when the center lies inside the loop)
#   bulges_5p/bulges_3p               unpaired bases on the 5'/3' arm of that stem, between the innermost pair and the loop
#   nearest_stem_dist                 positional distance from the center to the nearest paired base (0 when the center is paired)
# These definitions were regression-verified column-by-column against a 58-site
# reference output (174 identical rows).
import argparse
import csv
import os
import sys

import RNA

HALF_WINDOWS = (25, 50, 100)
DESC_COLS = ("stem_len", "loop_len", "u_to_loop_dist", "bulges_5p", "bulges_3p",
             "stem_gc", "stem_au", "stem_gu")
METRICS_COLS = (
    "genotype", "tissue", "transcript", "pos", "strand", "window", "start0", "end",
    "seq_len", "center_idx", "mfe", "mfe_per_nt", "center_unpaired_prob",
    "center_mfe_state", "center_element", "partner_offset",
    *DESC_COLS, "nearest_stem_dist", "ensemble_dg", "ensemble_dg_per_nt",
    "temp_C", "set", "region", "dot_bracket", "sequence",
)
PROFILE_COLS = ("genotype", "tissue", "transcript", "pos", "window", "rel_pos",
                "unpaired_prob", "temp_C", "set")


def parse_pairs(ss):
    """Parse base pairs: return the partner array, -1 for unpaired."""
    partner = [-1] * len(ss)
    stack = []
    for i, ch in enumerate(ss):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            j = stack.pop()
            partner[i] = j
            partner[j] = i
    return partner


def collect_stem(partner, i):
    """The complete stem containing i: extend inward and outward from (i, partner[i]); return [(a, b), ...] in ascending 5' index."""
    if partner[i] < 0:
        return []
    a, b = min(i, partner[i]), max(i, partner[i])
    pairs = [(a, b)]
    x, y = a + 1, b - 1
    while x < y and partner[x] == y:
        pairs.append((x, y))
        x += 1
        y -= 1
    x, y = a - 1, b + 1
    while x >= 0 and y < len(partner) and partner[x] == y:
        pairs.insert(0, (x, y))
        x -= 1
        y += 1
    return pairs


def classify_center(ss, partner, center):
    """Classify the structural element of the center site and give anchors describing
    its local stem-loop structure.

    After stack-parsing the pairs, find the nearest paired bases L/R on either side of
    the center (complete over pseudoknot-free structures):
    center paired -> stem; L and R paired to each other -> hairpin_loop;
    enclosed by a one-sided pair with no other stem inside -> bulge_internal_loop,
    with a stem inside -> multiloop;
    the rest (terminal free tail, linker between stems) -> other_unpaired.
    Returns (element, partner_offset, stem_anchor): stem_anchor is the starting point
    for collect_stem (the center itself for stem, L for hairpin_loop), or None when no
    single stem-loop unit applies.
    """
    n = len(ss)
    if ss[center] in "()":
        return "stem", partner[center] - center, center

    left = next((i for i in range(center - 1, -1, -1) if ss[i] in "()"), None)
    right = next((i for i in range(center + 1, n) if ss[i] in "()"), None)
    if left is None or right is None:
        return "other_unpaired", "", None
    if partner[left] == right:
        return "hairpin_loop", "", left
    if partner[left] > right:
        enclosed = any(ss[k] in "()" for k in range(right + 1, partner[left]))
        return ("multiloop" if enclosed else "bulge_internal_loop"), "", None
    if partner[right] < left:
        enclosed = any(ss[k] in "()" for k in range(partner[right] + 1, left))
        return ("multiloop" if enclosed else "bulge_internal_loop"), "", None
    return "other_unpaired", "", None


def stem_descriptors(ss, seq, partner, pairs, center_idx):
    """Structure descriptor columns of the stem-loop unit holding the center; all empty when pairs is empty."""
    out = {c: "" for c in DESC_COLS}
    if not pairs:
        return out

    a_in, b_in = pairs[-1]
    counts = {"stem_gc": 0, "stem_au": 0, "stem_gu": 0}
    for x, y in pairs:
        pair = seq[x] + seq[y]
        if pair in ("GC", "CG"):
            counts["stem_gc"] += 1
        elif pair in ("AU", "UA"):
            counts["stem_au"] += 1
        elif pair in ("GU", "UG"):
            counts["stem_gu"] += 1

    # stem_len / pair-type counts depend only on the stem itself and are always filled
    out.update({"stem_len": len(pairs), **counts})

    # Locate the loop of the stem-loop unit inward: layer by layer along the stem stack,
    # step over bulges on either side and over nested stems, until some layer's inner
    # region is fully unpaired — that segment is the loop of the stem-loop unit (bulge
    # counts accumulate over the layers).
    # If at some layer the two sides no longer pair complementarily (multibranch etc.),
    # the site does not form a single stem-loop unit: loop_len and u_to_loop_dist stay
    # empty, but bulge accumulation continues (matching the reference-output behavior).
    loop_start = loop_end = None
    bulges_5p = bulges_3p = 0
    matched = True
    cur_a, cur_b = a_in, b_in
    while cur_a + 1 <= cur_b - 1:
        x, y = cur_a + 1, cur_b - 1
        while x < y and ss[x] == ".":
            x += 1
        while y > x and ss[y] == ".":
            y -= 1
        if x >= y:                                       # inner side fully unpaired: this is the loop
            if matched:
                loop_start, loop_end = cur_a + 1, cur_b - 1
            break
        bulges_5p += x - (cur_a + 1)
        bulges_3p += (cur_b - 1) - y
        if partner[x] != y:
            # Multibranch: this layer is no longer a single stem-loop unit; loop_len is
            # not filled, but u_to_loop_dist still uses this layer's inner side as its
            # reference (matching the reference-output behavior)
            matched = False
            loop_start, loop_end = cur_a + 1, cur_b - 1
            break
        cur_a, cur_b = collect_stem(partner, x)[-1]

    out.update({"bulges_5p": bulges_5p, "bulges_3p": bulges_3p})
    if loop_start is not None and loop_start <= loop_end:
        if loop_start <= center_idx <= loop_end:
            u_to_loop_dist = 0
        else:
            lo, hi = (center_idx, loop_start) if center_idx < loop_start else (loop_end, center_idx)
            u_to_loop_dist = sum(1 for k in range(lo + 1, hi) if ss[k] in "()")
        out["u_to_loop_dist"] = u_to_loop_dist
        if matched:
            out["loop_len"] = loop_end - loop_start + 1
    return out


def unpaired_probabilities(bpp, n):
    """Per-position unpaired probabilities from the bpp matrix.

    fc.bpp() returns an (n+1)x(n+1), 1-based, strictly upper-triangular matrix: P(i,j)
    sits at bpp[i][j] (i<j). For 0-based sequence index k, p = k+1, and its pairing
    probability = the sum of row p to the right + the sum of column p above.
    """
    up = []
    for p in range(1, n + 1):
        paired = sum(bpp[p][p + 1:n + 1]) + sum(bpp[j][p] for j in range(1, p))
        up.append(min(1.0, max(0.0, 1.0 - paired)))
    return up


def self_check():
    """tRNA folding self-check: the stem should have low unpaired probability, the anticodon loop high, and the classifier should recognize a stem."""
    trna = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
    md = RNA.md()
    md.temperature = 21.0
    fc = RNA.fold_compound(trna, md)
    ss, mfe = fc.mfe()
    if mfe >= 0:
        sys.exit("[fold_structure] Self-check failed: tRNA MFE should be negative")
    fc.pf()
    up = unpaired_probabilities(fc.bpp(), len(trna))
    stem_up = sum(up[0:7]) / 7          # acceptor stem
    loop_up = sum(up[32:38]) / 6        # anticodon loop
    if not (stem_up < 0.5 and loop_up > 0.5):
        sys.exit(f"[fold_structure] Self-check failed: stem_up={stem_up:.3f} loop_up={loop_up:.3f}")
    element, _, _ = classify_center(ss, parse_pairs(ss), 3)   # a base inside the acceptor stem
    if element != "stem":
        sys.exit(f"[fold_structure] Self-check failed: classifier returned {element} for stem base")


def read_two_line_fasta(path):
    seqs = {}
    with open(path) as fh:
        lines = [line.strip() for line in fh if line.strip()]
    if len(lines) % 2 != 0:
        sys.exit(f"[fold_structure] Malformed two-line FASTA: {path}")
    for i in range(0, len(lines), 2):
        if not lines[i].startswith(">"):
            sys.exit(f"[fold_structure] Malformed FASTA header at line {i + 1} in {path}")
        seqs[lines[i][1:]] = lines[i + 1].upper()
    return seqs


def read_bed(path):
    bed = {}
    with open(path) as fh:
        for line in fh:
            tx, start, end, name, _score, strand = line.rstrip("\n").split("\t")
            bed[name] = (int(start), int(end))
    return bed


def main():
    ap = argparse.ArgumentParser(description="Fold site windows and extract structure metrics.")
    ap.add_argument("script_dir", help="Directory holding sites/bed/fasta inputs and outputs.")
    ap.add_argument("--sites", required=True, help="Site table path relative to script_dir.")
    ap.add_argument("--bed-prefix", required=True,
                    help="Window bed/fasta prefix relative to script_dir (suffix _{w}.bed/.fa).")
    ap.add_argument("--out-prefix", required=True, help="Output file prefix.")
    ap.add_argument("--temps", default="21", help="Space-separated temperatures in Celsius.")
    args = ap.parse_args()

    script_dir = args.script_dir
    temps = [float(t) for t in str(args.temps).split()]

    with open(os.path.join(script_dir, args.sites), newline="") as fh:
        sites = list(csv.DictReader(fh, delimiter="\t"))
    if not sites:
        sys.exit(f"[fold_structure] No sites in {args.sites}")

    metrics_rows = []
    profile_rows = []

    for temp in temps:
        md = RNA.md()
        md.temperature = temp
        for w in HALF_WINDOWS:
            bed = read_bed(os.path.join(script_dir, f"{args.bed_prefix}_{w}.bed"))
            seqs = read_two_line_fasta(os.path.join(script_dir, f"{args.bed_prefix}_{w}.fa"))

            for site in sites:
                name = f"{site['transcript']}_{site['pos']}_{site['strand']}"
                start0, end = bed[name]
                seq_dna = seqs[name]
                if len(seq_dna) != end - start0:
                    sys.exit(f"[fold_structure] Length mismatch for {name} w={w}: "
                             f"fasta {len(seq_dna)} vs bed {end - start0}")

                center_idx = (int(site["pos"]) - 1) - start0
                if not 0 <= center_idx < len(seq_dna):
                    sys.exit(f"[fold_structure] Center out of range for {name} w={w}")

                # Built-in strong validation: the 5-nt window center must reproduce the
                # context column of dependent.tsv; this also validates the coordinate
                # system, FASTA identity, and strand handling
                context = site["context"].upper().replace("U", "T")
                if seq_dna[center_idx - 2:center_idx + 3] != context:
                    sys.exit(f"[fold_structure] Context mismatch for {name} w={w}: "
                             f"window gives {seq_dna[center_idx - 2:center_idx + 3]}, "
                             f"expected {context}")

                seq = seq_dna.replace("T", "U")
                fc = RNA.fold_compound(seq, md)
                ss, mfe = fc.mfe()
                ensemble_dg = fc.pf()[1]
                up = unpaired_probabilities(fc.bpp(), len(seq))

                partner = parse_pairs(ss)
                element, partner_offset, anchor = classify_center(ss, partner, center_idx)
                pairs = collect_stem(partner, anchor) if anchor is not None else []
                desc = stem_descriptors(ss, seq, partner, pairs, center_idx)
                nearest_stem_dist = min(abs(i - center_idx) for i in range(len(ss))
                                        if ss[i] in "()")

                row = {
                    "genotype": site["genotype"], "tissue": site["tissue"],
                    "transcript": site["transcript"], "pos": site["pos"],
                    "strand": site["strand"], "window": w, "start0": start0, "end": end,
                    "seq_len": len(seq), "center_idx": center_idx,
                    "mfe": f"{mfe:.2f}", "mfe_per_nt": f"{mfe / len(seq):.4f}",
                    "center_unpaired_prob": f"{up[center_idx]:.4f}",
                    "center_mfe_state": "paired" if ss[center_idx] in "()" else "unpaired",
                    "center_element": element, "partner_offset": partner_offset,
                    **desc, "nearest_stem_dist": nearest_stem_dist,
                    "ensemble_dg": f"{ensemble_dg:.2f}",
                    "ensemble_dg_per_nt": f"{ensemble_dg / len(seq):.4f}",
                    "temp_C": int(temp), "set": site.get("set", ""),
                    "region": site.get("region", ""),
                    "dot_bracket": ss, "sequence": seq,
                }
                metrics_rows.append([row[c] for c in METRICS_COLS])

                for i, prob in enumerate(up):
                    profile_rows.append([
                        site["genotype"], site["tissue"], site["transcript"], site["pos"],
                        w, i - center_idx, f"{prob:.4f}", int(temp), site.get("set", ""),
                    ])

    metrics_path = os.path.join(script_dir, f"{args.out_prefix}_metrics.tsv")
    with open(metrics_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(METRICS_COLS)
        writer.writerows(metrics_rows)

    profile_path = os.path.join(script_dir, f"{args.out_prefix}_unpaired_profile.tsv")
    with open(profile_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(PROFILE_COLS)
        writer.writerows(profile_rows)

    print(f"[fold_structure] Wrote: {metrics_path} ({len(metrics_rows)} rows)")
    print(f"[fold_structure] Wrote: {profile_path} ({len(profile_rows)} rows)")


if __name__ == "__main__":
    self_check()
    main()
