#!/usr/bin/env python3
# build_case_sites.py — Build the case-site table for the trua1a-dependent
# secondary-structure comparison.
# Structurally identical to the case section of the sibling site_info/build_site_sets.py:
# the dependent rows of the three tissues are deduplicated into a transcript|pos|strand
# union, merging per-tissue information (tissues / n_tissues / max WT_Ratio, etc.);
# region comes from the union of the WT-positive site annotations of the three tissues
# (a transcript attribute, consistent across tissues); every case site must carry a
# complete annotation or the script fails.
# Only the case set is built: the controlPsi/controlUnmod sets of the sibling script
# serve case-control comparisons, which this figure set (union + per-tissue profiles)
# does not consume.
# Output columns are position-aligned: columns 3/4/5 are always transcript/pos/strand,
# compatible with the awk window-extraction logic of secondary_structure data_prep.sh.
import csv
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", "..", "..", "..", ".."))
DEPENDENT_DIR = f"{ROOT_DIR}/polyARNA/STAT/dependent_betabinomial/mRNA"
ANNO_DIR = f"{ROOT_DIR}/Figures/wildtype_landscape/mrna_region_distribution/site_analyse"

GENOTYPE = "trua1a"
TISSUES = ("Root", "leaf", "Flower")

OUT_CASE = f"{SCRIPT_DIR}/trua1a_dependent_sites.tsv"
OUT_REPORT = f"{SCRIPT_DIR}/site_sets_report.tsv"


def fail(msg):
    sys.exit(f"[build_case_sites] {msg}")


def key_of(transcript, pos, strand):
    return f"{transcript}|{pos}|{strand}"


def main():
    report = []  # (section, set, key, value)

    # ---------- case: unique physical trua1a-dependent sites ----------
    case = {}  # key -> dict
    n_rows = 0
    for tissue in TISSUES:
        path = f"{DEPENDENT_DIR}/{tissue}/{GENOTYPE}/{GENOTYPE}-{tissue}.dependent.tsv"
        if not os.path.isfile(path):
            fail(f"Missing input file: {path}")
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                n_rows += 1
                k = key_of(row["chr"], row["pos"], row["strand"])
                if k in case and case[k]["context"] != row["context"]:
                    fail(f"Case context mismatch across tissues: {k}")
                e = case.setdefault(k, {
                    "transcript": row["chr"], "pos": int(row["pos"]),
                    "strand": row["strand"], "context": row["context"],
                    "tissues": [], "WT_Ratio": 0.0, "Mut_Ratio": 0.0, "Delta": 0.0,
                })
                e["tissues"].append(tissue)
                e["WT_Ratio"] = max(e["WT_Ratio"], float(row["WT_Ratio"]))
                e["Mut_Ratio"] = max(e["Mut_Ratio"], float(row["Mut_Ratio"]))
                e["Delta"] = max(e["Delta"], float(row["Delta_Ratio(WT-Mut)"]))
    if not case:
        fail("No dependent sites found")
    report.append(("count", "case", "n_rows", n_rows))
    report.append(("count", "case", "n_unique_sites", len(case)))
    for n in (1, 2, 3):
        report.append(("count", "case", f"dependent_in_{n}_tissues",
                       sum(1 for e in case.values() if len(e["tissues"]) == n)))

    # ---------- WT-positive annotations (source of region) ----------
    region_by_key = {}   # key -> region (a transcript attribute, consistent across tissues)
    psi_tissues = {}     # key -> [tissue] (WT-positive union over the three tissues; counting only)
    for tissue in TISSUES:
        path = f"{ANNO_DIR}/WT-{tissue}_sites.positive.annoated.tsv"
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                k = key_of(row["chr"], row["pos"], row["strand"])
                if k in region_by_key and region_by_key[k] != row["annoation"]:
                    fail(f"Region mismatch across tissues: {k}")
                region_by_key[k] = row["annoation"]
                psi_tissues.setdefault(k, []).append(tissue)
    report.append(("count", "psi_union", "n_sites", len(psi_tissues)))

    # Every case site must lie inside the Ψ union (region annotation complete)
    for k, e in case.items():
        if k not in region_by_key:
            fail(f"Case site {k} lacks region annotation in positive.annoated.tsv")
        e["region"] = region_by_key[k]

    # ---------- Output (columns 3/4/5 are always transcript/pos/strand) ----------
    with open(OUT_CASE, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "tissues", "transcript", "pos", "strand", "context",
                    "region", "n_tissues", "WT_Ratio_max", "Mut_Ratio_max", "Delta_max"])
        for e in sorted(case.values(), key=lambda e: (e["transcript"], e["pos"], e["strand"])):
            w.writerow(["case", ",".join(e["tissues"]), e["transcript"], e["pos"],
                        e["strand"], e["context"], e["region"], len(e["tissues"]),
                        f"{e['WT_Ratio']:.6f}", f"{e['Mut_Ratio']:.6f}",
                        f"{e['Delta']:.6f}"])

    with open(OUT_REPORT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "set", "key", "value"])
        for row in report:
            w.writerow(row)

    print(f"[build_case_sites] case rows={n_rows} unique_sites={len(case)} "
          f"psi_union={len(psi_tissues)}")
    print(f"[build_case_sites] Wrote: {OUT_CASE}")
    print(f"[build_case_sites] Wrote: {OUT_REPORT}")


if __name__ == "__main__":
    main()
