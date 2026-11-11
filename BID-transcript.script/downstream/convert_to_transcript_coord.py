#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert genomic-coordinate BID sites to transcript-relative coordinates.

Background
----------
Genomic-coordinate-based pipelines identify sites using (chromosome, position, strand). 
This script maps these sites to the transcript-specific coordinate system used in the 
cDNA-alignment pipeline (where 'chr' is a transcript ID and 'pos' is a 1-based offset within the cDNA).

Modes
-----
--mode tx-col (default)
    Uses a transcript ID column within the input file to locate the site's transcript, 
    mapping coordinates based on that specific transcript's exon structure.

--mode genome
    Ignores the transcript ID column. Uses only genomic coordinates (chr, pos, strand) 
    to look up the site against a set of 'representative' transcripts (defined by --cdna .fai).
    Sites mapping to multiple representative transcripts are marked as 'ambiguous'.

Usage
-----
    python3 convert_to_transcript_coord.py --mode genome \
        --input sites.genome.tsv \
        --gtf <REPO_ROOT>/Araport11_fixed.gtf \
        --cdna <PATH_TO_REPRESENTATIVE_CDNA_COLLECTION> \
        --output sites.transcript.tsv

Output Columns
--------------
tx-col mode adds: `tx_pos`, `tx_pos_type`, `tx_base`, `tx_base_match`
genome mode adds: `rep_tx_id`, `rep_tx_pos`, `rep_tx_base`, `rep_base_match`,
                  `rep_assign`, `rep_n_candidates`, `rep_candidates`, `rep_spanning`
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict

BIN_SIZE = 10000


def normalize_chrom(value: str) -> str:
    """Standardize chromosome naming to Araport11_fixed.gtf format (Chr1..Chr5, ChrC, ChrM)."""
    text = str(value).strip()
    if len(text) == 1 and text in "12345":
        return "Chr" + text
    low = text.lower()
    if low in ("m", "mt", "chrm", "chrmt"):
        return "ChrM"
    if low in ("c", "pt", "chrc", "chrpt"):
        return "ChrC"
    return text


class Transcript:
    """Exon structure model for a single transcript."""

    __slots__ = ("tx_id", "chrom", "strand", "exons", "cdss", "length")

    def __init__(self, tx_id: str, chrom: str, strand: str):
        self.tx_id = tx_id
        self.chrom = chrom
        self.strand = strand
        self.exons: list[tuple[int, int]] = []
        self.cdss: list[tuple[int, int]] = []
        self.length = 0

    def build(self) -> None:
        self.exons.sort()
        self.cdss.sort()
        self.length = sum(end - start + 1 for start, end in self.exons)

    def in_cds(self, pos: int) -> bool:
        return any(start <= pos <= end for start, end in self.cdss)

    def gpos_to_tpos(self, pos: int) -> tuple[int, str]:
        """Convert genomic coordinate to transcript coordinate."""
        exons = self.exons
        if pos < exons[0][0]:
            dist = exons[0][0] - pos
            return dist, ("upstream" if self.strand == "+" else "downstream")
        if pos > exons[-1][1]:
            dist = pos - exons[-1][1]
            return dist, ("downstream" if self.strand == "+" else "upstream")

        tpos = 0
        for idx, (start, end) in enumerate(exons):
            if start <= pos:
                if end <= pos:
                    tpos += end - start + 1
                else:
                    tpos += pos - start + 1
            else:
                if exons[idx - 1][1] < pos:
                    if self.strand == "+":
                        return pos - exons[idx - 1][1], f"intron_{idx}"
                    return start - pos, f"intron_{len(exons) - idx}"
                break
        if self.strand == "-":
            tpos = self.length - tpos + 1
        return tpos, "exon"


class ExonIndex:
    """Exon range index supporting (chrom, pos, strand) lookup."""

    def __init__(self):
        self.bins: dict[tuple[str, int], list] = defaultdict(list)
        self.transcripts: dict[str, Transcript] = {}

    def add(self, tx: Transcript) -> None:
        self.transcripts[tx.tx_id] = tx
        order = tx.exons if tx.strand == "+" else tx.exons[::-1]
        offset = 0
        for start, end in order:
            record = (start, end, tx.tx_id, offset, tx.strand)
            for bucket in range(start // BIN_SIZE, end // BIN_SIZE + 1):
                self.bins[(tx.chrom, bucket)].append(record)
            offset += end - start + 1

    def query(self, chrom: str, pos: int, strand: str) -> list[tuple[str, int]]:
        hits = []
        for start, end, tx_id, offset, tx_strand in self.bins.get((chrom, pos // BIN_SIZE), ()):
            if tx_strand != strand or not (start <= pos <= end):
                continue
            if strand == "+":
                tx_pos = offset + (pos - start) + 1
            else:
                tx_pos = offset + (end - pos) + 1
            hits.append((tx_id, tx_pos))
        return hits

    def spanning(self, chrom: str, pos: int, strand: str) -> list[tuple[str, str]]:
        """Return mRNA range covering the position for representative transcripts."""
        found = []
        for tx in self.transcripts.values():
            if tx.chrom != chrom or tx.strand != strand:
                continue
            if tx.exons and tx.exons[0][0] <= pos <= tx.exons[-1][1]:
                tpos, ptype = tx.gpos_to_tpos(pos)
                found.append((tx.tx_id, ptype))
        return sorted(found)


def parse_gtf(path: str, keep: set[str] | None = None) -> ExonIndex:
    """Parse exon/CDS lines from GTF."""
    index = ExonIndex()
    pending: dict[str, Transcript] = {}
    cds_raw: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("exon", "CDS"):
                continue
            attrs = fields[8]
            key = 'transcript_id "'
            idx = attrs.find(key)
            if idx < 0:
                continue
            tx_id = attrs[idx + len(key):attrs.find('"', idx + len(key))]
            if keep is not None and tx_id not in keep:
                continue
            start, end = int(fields[3]), int(fields[4])
            if feature == "exon":
                tx = pending.get(tx_id)
                if tx is None:
                    tx = pending[tx_id] = Transcript(tx_id, fields[0], fields[6])
                tx.exons.append((start, end))
            else:
                cds_raw[tx_id].append((start, end))
    for tx_id, tx in pending.items():
        tx.cdss = cds_raw.get(tx_id, [])
        tx.build()
        index.add(tx)
    return index


def load_fai_names(path: str) -> set[str]:
    return {line.split("\t")[0] for line in open(path)}


def load_cdna(path: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def run_tx_col(args, index, cdna, header, rows) -> tuple[list[str], list[list[str]], dict]:
    tx_idx = header.index(args.tx_col)
    pos_idx = header.index(args.pos_col)
    tx_base_idx = header.index("transcript_base") if "transcript_base" in header else -1

    out_header = header + ["tx_pos", "tx_pos_type", "tx_base", "tx_base_match"]
    stats = defaultdict(int)
    out_rows = []
    for row in rows:
        tx_id, gpos = row[tx_idx], int(row[pos_idx])
        model = index.transcripts.get(tx_id)
        if model is None:
            tpos, ptype, base, match = 0, "NA", "", ""
        else:
            tpos, ptype = model.gpos_to_tpos(gpos)
            seq = cdna.get(tx_id)
            base = seq[tpos - 1] if seq and ptype == "exon" and 0 < tpos <= len(seq) else ""
            if tx_base_idx >= 0 and base:
                match = "TRUE" if base == row[tx_base_idx].upper() else "FALSE"
            else:
                match = ""
        stats[ptype if not ptype.startswith("intron") else "intron"] += 1
        if match:
            stats["base_" + match.lower()] += 1
        out_rows.append(row + [str(tpos), ptype, base, match])
    return out_header, out_rows, stats


def run_genome(args, index, cdna, header, rows) -> tuple[list[str], list[list[str]], dict]:
    chrom_idx = header.index(args.chrom_col)
    pos_idx = header.index(args.pos_col)
    strand_idx = header.index(args.strand_col) if args.strand_col in header else -1
    tx_base_idx = header.index("transcript_base") if "transcript_base" in header else -1

    out_header = header + [
        "rep_tx_id", "rep_tx_pos", "rep_tx_base", "rep_base_match",
        "rep_assign", "rep_n_candidates", "rep_candidates", "rep_spanning",
    ]
    stats = defaultdict(int)
    out_rows = []
    for row in rows:
        chrom = normalize_chrom(row[chrom_idx])
        gpos = int(row[pos_idx])
        strand = row[strand_idx] if strand_idx >= 0 else "+"

        hits = index.query(chrom, gpos, strand)
        antisense = index.query(chrom, gpos, "-" if strand == "+" else "+") if not hits else []

        if hits:
            if len(hits) == 1:
                tx_id, tpos = hits[0]
                assign = "unique"
            else:
                def rank(item):
                    cand_id, _ = item
                    cand = index.transcripts[cand_id]
                    return (0 if cand.in_cds(gpos) else 1, -cand.length, cand_id)
                tx_id, tpos = min(hits, key=rank)
                assign = "ambiguous"
            seq = cdna.get(tx_id)
            base = seq[tpos - 1] if seq and 0 < tpos <= len(seq) else ""
            if tx_base_idx >= 0 and base:
                match = "TRUE" if base == row[tx_base_idx].upper() else "FALSE"
            else:
                match = ""
            candidates = ";".join(f"{c}:{p}" for c, p in sorted(hits))
            stats["base_" + match.lower()] += 1 if match else 0
        else:
            tx_id, tpos, base, match = "", "", "", ""
            assign = "antisense_only" if antisense else "none"
            candidates = ""

        spanning = ""
        if not hits:
            spanning = ";".join(f"{t}:{p}" for t, p in index.spanning(chrom, gpos, strand))

        stats[assign] += 1
        out_rows.append(row + [
            tx_id, str(tpos), base, match, assign,
            str(len(hits)), candidates, spanning,
        ])
    return out_header, out_rows, stats


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Convert genomic-coordinate BID sites to transcript coordinates."
    )
    parser.add_argument("--input", required=True, help="Input TSV file with genomic coordinates.")
    parser.add_argument("--gtf", required=True, help="Araport11 GTF file.")
    parser.add_argument(
        "--cdna", default="",
        help="cDNA FASTA: used for self-check in tx-col mode; .fai used for representative transcripts in genome mode.",
    )
    parser.add_argument("--output", required=True, help="Output TSV file.")
    parser.add_argument("--mode", choices=("tx-col", "genome"), default="tx-col")
    parser.add_argument("--tx-col", default="transcript_id", help="Column name for transcript ID.")
    parser.add_argument("--pos-col", default="pos", help="Column name for genomic position.")
    parser.add_argument("--chrom-col", default="chr", help="Column name for chromosome.")
    parser.add_argument("--strand-col", default="strand", help="Column name for strand.")
    args = parser.parse_args()

    cdna = load_cdna(args.cdna) if args.cdna else {}
    if cdna:
        print(f"Loaded cDNA sequences: {len(cdna)}", file=sys.stderr)

    if args.mode == "genome":
        if not args.cdna:
            print("Mode 'genome' requires --cdna to determine representative transcripts.", file=sys.stderr)
            return 2
        keep = load_fai_names(args.cdna + ".fai")
        index = parse_gtf(args.gtf, keep=keep)
    else:
        index = parse_gtf(args.gtf)

    with open(args.input) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in fh]

    if args.mode == "genome":
        out_header, out_rows, stats = run_genome(args, index, cdna, header, rows)
    else:
        out_header, out_rows, stats = run_tx_col(args, index, cdna, header, rows)

    with open(args.output, "w") as fh:
        fh.write("\t".join(out_header) + "\n")
        for row in out_rows:
            fh.write("\t".join(row) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
