#!/usr/bin/env python
"""提取 WT mRNA Fisher 阳性位点三种窗口（±2/±5/±10）的序列表。

Fisher 阳性表（polyARNA/STAT/fisher/mRNA/WT-*_sites.positive.tsv）不含 context
列，因此 ±2/±5/±10 全部按转录本 1-based 拼接坐标从 Araport11 代表转录本 cDNA
（FASTA）切片：位点坐标为转录本 sense 方向，cDNA 同向，无需反向互补。
清洗规则与 dependent motif 管线一致：大写、U->T、只保留不含 N 的完整窗口序列。

每个窗口输出 <out_root>/+-<W>/mrna_motif_sequences.tsv（Tissue/sequence）
与 <out_root>/+-<W>/mrna_motif_center_base_stats.tsv（中心碱基分布 QC）。

Usage:
  python extract_window_sequences.py <fisher_mrna_dir> <fasta> <out_root>
"""

import csv
import sys
from pathlib import Path

WINDOWS = (2, 5, 10)
TISSUES = (("Root", "WT-Root"), ("Leaf", "WT-leaf"), ("Flower", "WT-Flower"))


def read_fasta(path: Path) -> dict:
    """读入 FASTA（id 取 header 首个 token），返回 {id: sequence}。"""
    sequences = {}
    name = None
    chunks = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            elif line:
                chunks.append(line)
    if name is not None:
        sequences[name] = "".join(chunks).upper()
    return sequences


def clean(seq: str) -> str:
    return seq.upper().replace("U", "T")


def main() -> None:
    fisher_dir = Path(sys.argv[1])
    fasta_path = Path(sys.argv[2])
    out_root = Path(sys.argv[3])

    fasta = read_fasta(fasta_path)
    print(f"[extract] FASTA transcripts: {len(fasta)}")

    for window in WINDOWS:
        print(f"[extract] Window ±{window}")
        window_dir = out_root / f"+-{window}"
        window_dir.mkdir(parents=True, exist_ok=True)
        stats_rows = []

        seq_file = window_dir / "mrna_motif_sequences.tsv"
        with seq_file.open("w", newline="") as seq_handle:
            writer = csv.writer(seq_handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Tissue", "sequence"])

            for label, prefix in TISSUES:
                input_file = fisher_dir / f"{prefix}_sites.positive.tsv"
                if not input_file.is_file():
                    raise SystemExit(f"Missing input: {input_file}")

                sequences = []
                n_input = n_drop_tx = n_drop_edge = n_drop_n = 0
                half = window

                with input_file.open(newline="") as handle:
                    reader = csv.DictReader(handle, delimiter="\t")
                    required = {"chr", "pos"}
                    if not required.issubset(reader.fieldnames or []):
                        raise SystemExit(f"Unexpected columns in {input_file}: {reader.fieldnames}")
                    for row in reader:
                        n_input += 1
                        tx = row["chr"]
                        pos = int(row["pos"])
                        transcript = fasta.get(tx)
                        if transcript is None:
                            n_drop_tx += 1
                            continue
                        start0 = pos - 1 - half
                        end0 = pos + half
                        if start0 < 0 or end0 > len(transcript):
                            n_drop_edge += 1
                            continue
                        seq = clean(transcript[start0:end0])
                        if set(seq) - set("ACGT"):
                            n_drop_n += 1
                            continue
                        sequences.append(seq)

                writer.writerows((label, seq) for seq in sequences)

                counts = {base: 0 for base in "ACGT"}
                for seq in sequences:
                    counts[seq[half]] += 1
                stats_rows.append(
                    [label, counts["A"], counts["C"], counts["G"], counts["T"], sum(counts.values())]
                )
                print(
                    f"[extract]   {label}: input={n_input} kept={len(sequences)} "
                    f"drop_tx_absent={n_drop_tx} drop_edge={n_drop_edge} drop_invalid={n_drop_n}"
                )

        stats_file = window_dir / "mrna_motif_center_base_stats.tsv"
        with stats_file.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["Tissue", "A", "C", "G", "T", "Total"])
            writer.writerows(stats_rows)


if __name__ == "__main__":
    main()
