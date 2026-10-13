#!/usr/bin/env python3
# build_site_sets.py — 为 trua1a dependent 精细化二级结构对比构建三组位点
# 职责：
#   1) case：trua1a dependent 三组织全部位点行 → unique 物理位点（transcript|pos|strand
#      去重），合并各组织 dependent 信息（tissues/n_tissues/最大 WT_Ratio 等）；
#   2) controlPsi（TRUA1A 酶催化特异性对照）：三组织 WT 阳性 Ψ 位点并集排除
#      case 位点；context 由参考 FASTA 随机访问提取并校验中心碱基为 T；
#   3) controlUnmod（未修饰背景对照）：按 Pseudo_callsites.R 同款规则从 filter_T.tsv
#      （参考碱基 T 且 depth≥10）重建各样品可测池——全部重复 merge(all=TRUE) 缺失补 0，
#      组均 depth≥10——不做
#      p.adj 信号过滤；取全部样品池交集后要求每个样品 Treat_avg_ratio < 0.02。
#      （2026-09-15 修正：不可从 Fisher.tsv 选未修饰位点——Fisher.tsv 写出前已过滤
#      pooled p.adj<0.01，全是 Treat 显著大于 Mock 的信号位点，并非未修饰宇宙。）
# 输出三表列位对齐：第 3/4/5 列恒为 transcript/pos/strand，兼容 secondary_structure
# data_prep.sh 式 awk 窗口提取逻辑；第 2 列 tissues 语义：case/controlPsi 为鉴定到该位点
# 的组织列表，controlUnmod 恒为 Root,leaf,Flower（45 样品交集的必然结果）。
# 另出 site_sets_report.tsv 记录各步计数；三集合两两不重叠，重叠即报错。
import csv
import os
import sys

import polars as pl

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "..", "..", "..", "..", "..", ".."))
DEPENDENT_DIR = f"{ROOT_DIR}/polyARNA/STAT/dependent_methylKit/mRNA"
FILTER_T_DIR = f"{ROOT_DIR}/polyARNA/STAT/adjustGap/mRNA"
ANNO_DIR = f"{ROOT_DIR}/Figures/wildtype_landscape/mrna_region_distribution/site_analyse"
REF_FA = "/home/HL/reference/Ath/TAIR/Araport11_cdna_20220914_representative_gene_model"

GENOTYPE = "trua1a"
TISSUES = ("Root", "leaf", "Flower")
RATIO_CAP = 0.02     # controlUnmod：所有样品 Treat_avg_ratio 上限
MIN_AVG_DEPTH = 10   # Pseudo_callsites.R --min-avg-depth 默认阈值
# 分析口径固定为每样品 3 个生物学重复（WT 的第 4 重复来自不同批次、
# 批次效应过大，已于 2026-09-26 弃用）。WT 的 rep4 filter_T 文件仍留在
# adjustGap 目录里（上游不重跑），故此处必须显式排除，否则可测池交集
# 会按 4 重复口径计算。
MAX_REPS = 3
KEY4 = ["chr", "pos", "strand", "context"]

OUT_CASE = f"{SCRIPT_DIR}/trua1a_dependent_sites.tsv"
OUT_PSI = f"{SCRIPT_DIR}/control_nondependent_psi_sites.tsv"
OUT_UNMOD = f"{SCRIPT_DIR}/control_unmodified_sites.tsv"
OUT_REPORT = f"{SCRIPT_DIR}/site_sets_report.tsv"


def fail(msg):
    sys.exit(f"[build_site_sets] {msg}")


def key_of(transcript, pos, strand):
    return f"{transcript}|{pos}|{strand}"


class FaiReader:
    """基于 .fai 的随机访问 FASTA 读取器（处理固定行宽的换行偏移）。"""

    def __init__(self, fa_path):
        if not os.path.isfile(fa_path) or not os.path.isfile(fa_path + ".fai"):
            fail(f"Missing reference FASTA or .fai: {fa_path}")
        self.fh = open(fa_path)
        self.index = {}
        with open(fa_path + ".fai") as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                self.index[f[0]] = (int(f[1]), int(f[2]), int(f[3]), int(f[4]))

    def fetch(self, name, start0, end0):
        """0-based 半开区间切片。"""
        if name not in self.index:
            fail(f"Transcript {name} not in reference .fai")
        length, offset, line_bases, line_width = self.index[name]
        if start0 < 0 or end0 > length:
            fail(f"fetch out of range: {name}:{start0}-{end0} (len {length})")
        buf = []
        for i in range(start0, end0):
            self.fh.seek(offset + (i // line_bases) * line_width + i % line_bases)
            buf.append(self.fh.read(1))
        return "".join(buf).upper()


def sample_measurable_lf(sample, rep_files):
    """单样品可测池 LazyFrame：该样品全部 Mock/Treat 重复 full join（缺失补 0）→
    组均 depth/ratio → Mock/Treat 组均 depth >= MIN_AVG_DEPTH（Pseudo_callsites.R
    规则，无 p.adj 过滤）。重复数按样品实际持有的文件动态确定：突变体系 3 重复，
    WT 自 2608 批次补测后为 4 重复。"""
    lf = None
    for cond, rep, path in rep_files:
        tag = f"{cond}{rep}"
        sub = (pl.scan_csv(path, separator="\t",
                           schema_overrides={"chr": pl.String, "pos": pl.Int64,
                                             "strand": pl.String, "context": pl.String,
                                             "ratio": pl.Float64, "depth": pl.Float64})
               .select(KEY4 + [pl.col("ratio").alias(f"{tag}_ratio"),
                               pl.col("depth").alias(f"{tag}_depth")]))
        lf = sub if lf is None else lf.join(sub, on=KEY4, how="full", coalesce=True)
    value_cols = [c for c in lf.collect_schema().names() if c not in KEY4]
    mock_reps = sorted(rep for cond, rep, _ in rep_files if cond == "Mock")
    treat_reps = sorted(rep for cond, rep, _ in rep_files if cond == "Treat")
    mock_depths = [f"Mock{r}_depth" for r in mock_reps]
    treat_depths = [f"Treat{r}_depth" for r in treat_reps]
    treat_ratios = [f"Treat{r}_ratio" for r in treat_reps]
    return (lf.with_columns(pl.col(c).fill_null(0.0) for c in value_cols)
            .with_columns(Mock_avg_depth=pl.mean_horizontal(mock_depths),
                          Treat_avg_depth=pl.mean_horizontal(treat_depths),
                          Treat_avg_ratio=pl.mean_horizontal(treat_ratios))
            .filter((pl.col("Mock_avg_depth") >= MIN_AVG_DEPTH)
                    & (pl.col("Treat_avg_depth") >= MIN_AVG_DEPTH))
            .select(KEY4 + ["Treat_avg_ratio"]))


def main():
    report = []  # (section, set, key, value)

    # ---------- 1) case：trua1a dependent unique 物理位点 ----------
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

    # ---------- WT 阳性注释（region 来源 + Ψ 位点并集） ----------
    region_by_key = {}   # key -> region（转录本属性，跨组织一致）
    psi_tissues = {}     # key -> [tissue]（三组织 WT 阳性并集）
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

    # 病例必须全部落在 Ψ 并集内（region 注释齐全）
    for k, e in case.items():
        if k not in region_by_key:
            fail(f"Case site {k} lacks region annotation in positive.annoated.tsv")
        e["region"] = region_by_key[k]

    # ---------- 2) controlPsi：Ψ 并集排除 dependent，FASTA 补 context ----------
    fai = FaiReader(REF_FA)
    control_psi = {}
    for k, tissues in psi_tissues.items():
        if k in case:
            continue
        tx, pos, strand = k.split("|")
        pos = int(pos)
        context = fai.fetch(tx, pos - 3, pos + 2)  # pos 1-based → 中心 5-mer
        if len(context) != 5 or context[2] != "T":
            fail(f"controlPsi center base is not T: {k} gives {context}")
        control_psi[k] = {
            "transcript": tx, "pos": pos, "strand": strand, "context": context,
            "region": region_by_key[k], "tissues": tissues,
        }
    report.append(("count", "controlPsi", "n_sites", len(control_psi)))
    report.append(("excluded", "controlPsi", "case_sites_removed",
                   len(psi_tissues) - len(control_psi)))

    # ---------- 3) controlUnmod：filter_T 重建样品可测池 → 交集 + ratio 上限 ----------
    # 样品清单：Mock/Treat-<sample>.rep*.mRNA.filter_T.tsv，重复数按实际文件动态确定，
    # 但一律截断到 MAX_REPS（突变体系本就 3 个；WT 磁盘上有 rep1-4，须排除 rep4）。
    sample_files = {}
    for fname in sorted(os.listdir(FILTER_T_DIR)):
        if not fname.endswith(".mRNA.filter_T.tsv"):
            continue
        stem = fname[: -len(".mRNA.filter_T.tsv")]  # Mock-WT-Root.rep1
        cond_sample, rep = stem.rsplit(".rep", 1)
        if int(rep) > MAX_REPS:
            continue
        cond, sample = cond_sample.split("-", 1)
        sample_files.setdefault(sample, []).append(
            (cond, int(rep), f"{FILTER_T_DIR}/{fname}"))
    for sample, fl in sorted(sample_files.items()):
        mock_reps = sorted(r for c, r, _ in fl if c == "Mock")
        treat_reps = sorted(r for c, r, _ in fl if c == "Treat")
        if not mock_reps or mock_reps != treat_reps:
            fail(f"Sample {sample} has unmatched Mock/Treat filter_T files: "
                 f"Mock rep{mock_reps} vs Treat rep{treat_reps}")
        fl.sort(key=lambda x: (x[0], x[1]))
    n_samples = len(sample_files)
    report.append(("count", "controlUnmod", "n_samples", n_samples))

    acc = None  # 交集累积帧：KEY4 + ratio_max + ratio_sum + n_present
    for i, (sample, fl) in enumerate(sorted(sample_files.items()), start=1):
        df = sample_measurable_lf(sample, fl).collect(engine="streaming")
        report.append(("measurable", "controlUnmod", sample, df.height))
        if acc is None:
            acc = df.rename({"Treat_avg_ratio": "ratio_max"}).with_columns(
                ratio_sum=pl.col("ratio_max"), n_present=pl.lit(1))
        else:
            acc = (acc.join(df, on=KEY4, how="inner")
                   .with_columns(ratio_max=pl.max_horizontal("ratio_max", "Treat_avg_ratio"),
                                 ratio_sum=pl.col("ratio_sum") + pl.col("Treat_avg_ratio"),
                                 n_present=pl.col("n_present") + 1)
                   .drop("Treat_avg_ratio"))
        report.append(("intersection", "controlUnmod", sample, acc.height))
        print(f"[build_site_sets] ({i}/{n_samples}) {sample}: measurable={df.height} "
              f"intersection={acc.height}", flush=True)
        del df
    if acc is None or acc.height == 0:
        fail("Empty intersection across samples")
    if acc.filter(pl.col("n_present") != n_samples).height:  # inner join 已保证，防御
        fail("Intersection contains sites not present in all samples")
    n_intersect = acc.height
    # 转录本 5' 端边缘位点（pos<=2）context 为截断短串、中心非 T，无法做 5-mer/对称窗口，
    # 剔除并计数（参考碱基 T 由中心碱基校验同时确认）
    edge = acc.filter((pl.col("context").str.len_chars() != 5)
                      | (pl.col("context").str.slice(2, 1) != "T")).height
    acc = acc.filter((pl.col("context").str.len_chars() == 5)
                     & (pl.col("context").str.slice(2, 1) == "T"))
    report.append(("count", "controlUnmod", "n_intersect_sites", n_intersect))
    report.append(("excluded", "controlUnmod", "edge_or_bad_context", edge))
    unmod_df = (acc.filter(pl.col("ratio_max") < RATIO_CAP)
                .with_columns(Treat_ratio_mean=pl.col("ratio_sum") / n_samples)
                .select(KEY4 + ["ratio_max", "Treat_ratio_mean"]))
    n_unmod = unmod_df.height
    report.append(("count", "controlUnmod", f"n_ratio_cap_{RATIO_CAP}", n_unmod))
    if n_unmod == 0:
        fail("No unmodified sites survive the ratio cap")
    del acc

    # region 注释：三组织 WT 表达池（measurable-in-all-45 ⊆ 各 WT 池，理论上全覆盖）
    pool = pl.concat([pl.scan_csv(f"{ANNO_DIR}/WT-{t}.expressed_pool.annoated.tsv",
                                  separator="\t")
                      .select("chr", "pos", "strand", region="annoation")
                      for t in TISSUES]).unique(subset=["chr", "pos", "strand"],
                                                 keep="first")
    unmod_df = (unmod_df.join(pool.collect(engine="streaming"),
                              on=["chr", "pos", "strand"], how="left")
                .with_columns(pl.col("region").fill_null("unmapped")))
    n_no_region = unmod_df.filter(pl.col("region") == "unmapped").height
    report.append(("count", "controlUnmod", "n_region_annotated", n_unmod - n_no_region))
    report.append(("count", "controlUnmod", "n_region_missing", n_no_region))

    # ---------- 三集合两两不重叠 ----------
    unmod_keys = set(unmod_df.select(
        k=pl.concat_str("chr", "pos", "strand", separator="|")).to_series())
    overlap_cp = set(case) & set(control_psi)
    overlap_cu = set(case) & unmod_keys
    overlap_pu = set(control_psi) & unmod_keys
    if overlap_cp or overlap_cu or overlap_pu:
        fail(f"Set overlap: case∩psi={len(overlap_cp)} case∩unmod={len(overlap_cu)} "
             f"psi∩unmod={len(overlap_pu)}")

    # ---------- 写出（第 3/4/5 列恒为 transcript/pos/strand） ----------
    def sort_key(e):
        return (e["transcript"], e["pos"], e["strand"])

    with open(OUT_CASE, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "tissues", "transcript", "pos", "strand", "context",
                    "region", "n_tissues", "WT_Ratio_max", "Mut_Ratio_max", "Delta_max"])
        for e in sorted(case.values(), key=sort_key):
            w.writerow(["case", ",".join(e["tissues"]), e["transcript"], e["pos"],
                        e["strand"], e["context"], e["region"], len(e["tissues"]),
                        f"{e['WT_Ratio']:.6f}", f"{e['Mut_Ratio']:.6f}",
                        f"{e['Delta']:.6f}"])

    with open(OUT_PSI, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["set", "tissues", "transcript", "pos", "strand", "context",
                    "region", "n_tissues"])
        for e in sorted(control_psi.values(), key=sort_key):
            w.writerow(["controlPsi", ",".join(e["tissues"]), e["transcript"], e["pos"],
                        e["strand"], e["context"], e["region"], len(e["tissues"])])

    (unmod_df.sort("chr", "pos", "strand")
     .with_columns(set=pl.lit("controlUnmod"), tissues=pl.lit(",".join(TISSUES)),
                   n_samples=pl.lit(n_samples))  # ratio 列保留全精度：四舍五入会把
     .select(pl.col("set"), pl.col("tissues"), pl.col("chr").alias("transcript"),  # 0.0199996
             pl.col("pos"), pl.col("strand"), pl.col("context"), pl.col("region"),  # 显示成 0.02
             pl.col("n_samples"), pl.col("ratio_max").alias("Treat_ratio_max"),
             pl.col("Treat_ratio_mean"))
     .write_csv(OUT_UNMOD, separator="\t"))

    with open(OUT_REPORT, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["section", "set", "key", "value"])
        w.writerows(report)

    print(f"[build_site_sets] case={len(case)} controlPsi={len(control_psi)} "
          f"controlUnmod={n_unmod} (intersection={n_intersect})")
    print(f"[build_site_sets] Wrote: {OUT_CASE}")
    print(f"[build_site_sets] Wrote: {OUT_PSI}")
    print(f"[build_site_sets] Wrote: {OUT_UNMOD}")
    print(f"[build_site_sets] Wrote: {OUT_REPORT}")


if __name__ == "__main__":
    main()
