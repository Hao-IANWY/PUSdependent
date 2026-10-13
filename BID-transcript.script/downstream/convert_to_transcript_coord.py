#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把基因组坐标的 BID-seq 位点表转换为转录组坐标。

背景
----
另一套流程在**基因组**上比对调用位点，输出的表用 染色体 + 基因组坐标 + 链向 定位；
本仓库的流程（BID-transcript.script / polyA_mapping.smk）把 reads 直接比对到**转录组**
cDNA 参考，`STAT/` 下所有表的 `chr` 列是转录本 ID、`pos` 是该转录本内的 1-based 坐标。
两套坐标无法直接比对，需要本脚本换算。

两种模式
--------
`--mode tx-col`（默认）
    用表内的转录本 ID 列指定位点所属转录本，再按该转录本的外显子结构换算。
    允许是非代表异构体（只要 GTF 里有它的外显子结构）。

`--mode genome`
    完全忽略表内的转录本 ID 列，只用 (chr, pos, strand) 基因组坐标，
    在本流程实际使用的**代表转录本**集合里反查包含该位点的转录本。
    代表转录本集合由 `--cdna` 的 .fai 定义（即流程比对用的那套 cDNA）。
    同一坐标落在多条代表转录本上时（基因重叠）取主选并标记 `ambiguous`。

实现不依赖外部 gppy 工具，直接在脚本内解析 GTF 完成映射。

用法
----
    python3 convert_to_transcript_coord.py --mode genome \
        --input  sites.genome.tsv \
        --gtf    <REPO_ROOT>/Araport11_fixed.gtf \
        --cdna   /home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model \
        --output sites.transcript.tsv

输出列
------
tx-col 模式追加：`tx_pos` / `tx_pos_type` / `tx_base` / `tx_base_match`
genome 模式追加：`rep_tx_id` / `rep_tx_pos` / `rep_tx_base` / `rep_base_match`
                 / `rep_assign` / `rep_n_candidates` / `rep_candidates`
"""

from __future__ import annotations

import argparse
import sys
from collections import defaultdict

BIN_SIZE = 10000


def normalize_chrom(value: str) -> str:
    """统一染色体命名到 Araport11_fixed.gtf 的口径（Chr1..Chr5 / ChrC / ChrM）。"""
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
    """一条转录本的外显子模型。"""

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
        """基因组坐标 -> (转录本坐标, 类型)。语义与 gppy gtf.py 的 gpos_to_tpos 一致。

        转录本坐标从 5' 端起算（负链即基因组高坐标端）。返回值：
          - `exon`：正常的外显子坐标，取值 1..转录本长度
          - `intron_N`：落在内含子内（N 从 5' 端起数），坐标是该内含子内的偏移
          - `upstream` / `downstream`：落在整个 mRNA 跨度之外，坐标是到边界的碱基数
        """
        exons = self.exons  # 基因组升序
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
                # pos 落在本外显子之前：若不在上一外显子内，则处于两者之间的内含子
                if exons[idx - 1][1] < pos:
                    if self.strand == "+":
                        return pos - exons[idx - 1][1], f"intron_{idx}"
                    return start - pos, f"intron_{len(exons) - idx}"
                break
        if self.strand == "-":
            tpos = self.length - tpos + 1
        return tpos, "exon"


class ExonIndex:
    """外显子区间索引，支持 (chrom, pos, strand) 反查转录本。"""

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
        """反查 mRNA 基因组跨度覆盖该位置的代表转录本，返回 (tx_id, 位置类型)。

        用于解释外显子反查落空的原因：位点常常落在代表转录本的内含子里。
        """
        found = []
        for tx in self.transcripts.values():
            if tx.chrom != chrom or tx.strand != strand:
                continue
            if tx.exons and tx.exons[0][0] <= pos <= tx.exons[-1][1]:
                tpos, ptype = tx.gpos_to_tpos(pos)
                found.append((tx.tx_id, ptype))
        return sorted(found)


def parse_gtf(path: str, keep: set[str] | None = None) -> ExonIndex:
    """解析 GTF 的 exon/CDS 行。keep 非空时只保留其中的转录本。"""
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
    """.fai 的序列名即转录本 ID。"""
    return {line.split("\t")[0] for line in open(path)}


def load_cdna(path: str) -> dict[str, str]:
    """读取 cDNA FASTA（仅用于自检碱基），返回 {transcript_id: sequence}。"""
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
    """模式一：按表内转录本 ID 列换算。"""
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
    """模式二：只用基因组坐标，在代表转录本里反查。"""
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
                # 主选规则：优先落在 CDS 内，其次取更长的转录本，最后按 ID 字典序
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

        # 说明性列：外显子反查落空时，列出 mRNA 跨度覆盖该位置的代表转录本及其类型
        # (通常显示 intron_N —— 该位点在代表异构体里是内含子, 即选择性外显子)
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
    parser.add_argument("--input", required=True, help="输入的基因组坐标位点表（TSV）")
    parser.add_argument("--gtf", required=True, help="Araport11 GTF")
    parser.add_argument(
        "--cdna", default="",
        help="cDNA FASTA：tx-col 模式用于生成自检列；genome 模式另用其 .fai 定义代表转录本集合",
    )
    parser.add_argument("--output", required=True, help="输出 TSV")
    parser.add_argument("--mode", choices=("tx-col", "genome"), default="tx-col")
    parser.add_argument("--tx-col", default="transcript_id", help="转录本 ID 列名")
    parser.add_argument("--pos-col", default="pos", help="基因组坐标列名")
    parser.add_argument("--chrom-col", default="chr", help="染色体列名")
    parser.add_argument("--strand-col", default="strand", help="链向列名")
    args = parser.parse_args()

    cdna = load_cdna(args.cdna) if args.cdna else {}
    if cdna:
        print(f"[convert] cDNA 序列数: {len(cdna)}", file=sys.stderr)

    if args.mode == "genome":
        if not args.cdna:
            print("[convert] genome 模式需要 --cdna 以确定代表转录本集合", file=sys.stderr)
            return 2
        keep = load_fai_names(args.cdna + ".fai")
        index = parse_gtf(args.gtf, keep=keep)
        print(
            f"[convert] 代表转录本 {len(keep)} 条，GTF 中命中 {len(index.transcripts)} 条",
            file=sys.stderr,
        )
    else:
        index = parse_gtf(args.gtf)
        print(f"[convert] GTF 转录本数: {len(index.transcripts)}", file=sys.stderr)

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

    print(f"[convert] 共 {len(rows)} 个位点", file=sys.stderr)
    for key in ("unique", "ambiguous", "antisense_only", "none", "exon", "intron", "NA"):
        if stats.get(key):
            print(f"[convert]   {key}: {stats[key]}", file=sys.stderr)
    if stats.get("base_true") or stats.get("base_false"):
        print(
            f"[convert] 碱基自检: 一致 {stats.get('base_true', 0)}, "
            f"不一致 {stats.get('base_false', 0)}",
            file=sys.stderr,
        )
    print(f"[convert] 输出: {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
