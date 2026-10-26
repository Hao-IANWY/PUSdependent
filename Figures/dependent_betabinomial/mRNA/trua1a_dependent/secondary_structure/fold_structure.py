#!/usr/bin/env python3
# fold_structure.py — RNAfold 折叠位点窗口序列并提取结构指标
#
# 职责：对 sites 表里的每个位点 × 每档窗口(±25/±50/±100) × 每个温度，用 ViennaRNA
#       计算 MFE 结构与配对概率，分类中心位点的结构元件，输出 metrics 与逐位
#       unpaired profile 两个 TSV。
#
# 用法：
#   fold_structure.py <script_dir> --sites <rel> --bed-prefix <rel> \
#                     --out-prefix <name> [--temps "21 37"]
#   输入 —— <script_dir>/<sites>：genotype/tissue/transcript/pos/strand/context[/set]
#           <script_dir>/<bed-prefix>_{w}.bed + _{w}.fa：bedtools ±w 窗口
#   输出 —— <script_dir>/<out-prefix>_metrics.tsv
#           <script_dir>/<out-prefix>_unpaired_profile.tsv
#
# 模型参数：ViennaRNA 默认（dangles=2、允许 GU 摆动配对、无伪结），温度由 --temps 给定
#           （本项目正文分析用 21 °C 体内近似）。
#
# 结构描述列（仅当中心归属单一茎环单元时非空）：
#   stem_len/stem_gc/stem_au/stem_gu  中心所在茎的配对数与配对类型计数
#   loop_len                          该茎围成的最内层环长（茎内层对之间的未配对段）
#   u_to_loop_dist                    中心到该环之间跨过的配对碱基数（中心在环内为 0）
#   bulges_5p/bulges_3p               该茎 5'/3' 链上、内层对到环之间的未配对碱基数
#   nearest_stem_dist                 中心到最近配对碱基的位置距离（中心配对时为 0）
# 上述定义经 2026-09-24 用 58 位点旧产出逐列回归验证（174 行全等）。
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
    """碱基配对解析：返回 partner 数组，未配对为 -1。"""
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
    """含 i 的完整茎：从 (i, partner[i]) 向内外两侧扩展，按 5' 索引升序返回 [(a, b), ...]。"""
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
    """分类中心位点的结构元件，并给出描述该位点局部茎环结构的锚点。

    栈解析配对后从中心向两侧找最近的配对碱基 L/R（伪结自由结构下情形完备）：
    中心已配对 -> stem；L/R 互为配对 -> hairpin_loop；
    被一侧配对包住、内侧无其他茎 -> bulge_internal_loop，有茎 -> multiloop；
    其余（端部自由尾、茎间 linker）-> other_unpaired。
    返回 (element, partner_offset, stem_anchor)：stem_anchor 为 collect_stem 的起点
    （stem 用中心本身，hairpin_loop 用 L），无单一茎环单元时为 None。
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
    """中心所在茎环单元的结构描述列；pairs 为空时全部留空。"""
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

    # stem_len / 配对类型计数只依赖茎本身，恒定填写
    out.update({"stem_len": len(pairs), **counts})

    # 向内定位该茎环单元的环：沿茎堆叠逐层跨过两侧 bulge 与被套住的茎，直到某层的
    # 内侧区域完全未配对——那一段就是该茎环单元的环（bulge 数按各层累积）。
    # 途中若某层两侧不再互补配对（多分支等），该位点不构成单一茎环单元：loop_len 与
    # u_to_loop_dist 留空，但 bulge 累积照常继续（与旧产出行为一致）。
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
        if x >= y:                                       # 内侧全未配对：即环
            if matched:
                loop_start, loop_end = cur_a + 1, cur_b - 1
            break
        bulges_5p += x - (cur_a + 1)
        bulges_3p += (cur_b - 1) - y
        if partner[x] != y:
            # 多分支：该层不再是单一茎环单元，loop_len 不再填写，但 u_to_loop_dist
            # 仍以该层内侧为参照给出（与旧产出行为一致）
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
    """由 bpp 矩阵算逐位 unpaired probability。

    fc.bpp() 返回 (n+1)x(n+1)、1-based、严格上三角矩阵：P(i,j) 在 bpp[i][j] (i<j)。
    0-based 序列下标 k 对应 p = k+1，其配对概率 = 行 p 右侧之和 + 列 p 上方之和。
    """
    up = []
    for p in range(1, n + 1):
        paired = sum(bpp[p][p + 1:n + 1]) + sum(bpp[j][p] for j in range(1, p))
        up.append(min(1.0, max(0.0, 1.0 - paired)))
    return up


def self_check():
    """tRNA 折叠自检：茎区 unpaired 概率应低、反密码子环应高，分类器应识别 stem。"""
    trna = "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUCUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCACCA"
    md = RNA.md()
    md.temperature = 21.0
    fc = RNA.fold_compound(trna, md)
    ss, mfe = fc.mfe()
    if mfe >= 0:
        sys.exit("[fold_structure] Self-check failed: tRNA MFE should be negative")
    fc.pf()
    up = unpaired_probabilities(fc.bpp(), len(trna))
    stem_up = sum(up[0:7]) / 7          # 受体茎
    loop_up = sum(up[32:38]) / 6        # 反密码子环
    if not (stem_up < 0.5 and loop_up > 0.5):
        sys.exit(f"[fold_structure] Self-check failed: stem_up={stem_up:.3f} loop_up={loop_up:.3f}")
    element, _, _ = classify_center(ss, parse_pairs(ss), 3)   # 受体茎内碱基
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

                # 内置强校验：窗口中心 5nt 必须复现 dependent.tsv 的 context 列，
                # 同时验证坐标系、FASTA 身份与链处理
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
