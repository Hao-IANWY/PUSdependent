# tRNA Sprinzl 结构注释表

## 目的

本目录保存 tRNA 位点坐标转换（native → canonical / Sprinzl 0–76 位）所依据的两张注释表：

- **表一** `plantRNA_database.tsv`：从 plantRNA 数据库下载的原始注释表（644 条基因组拷贝记录）。
- **表二** `plantRNA_database_unique.tsv`：对表一的序列取 unique 后得到的 220 条 unique tRNA 注释，与比对参考里的 220 条 tRNA 一一对应。

下游用法见 `Figures/dependent_methylKit/tRNA/build_genotype_dependent_tables.sh`（把 native pos 按区段重标定为 `fix_pos`）。

## 目录内容

```
tRNA_Sprinzl_annoation/
├── plantRNA_database.tsv            表一：plantRNA 数据库下载的注释（644 条基因组拷贝）
├── plantRNA_database_unique.tsv     表二：220 条 unique tRNA（每条一份结构区段注释）
├── add_sprinzl_pos.py               给位点表追加 Sprinzl_pos 列的脚本
└── README.md
```

## 表一：plantRNA_database.tsv

- **来源**：plantRNA 数据库下载的注释表。下载原件为 `.xls`，本 TSV 是**按原列名、原取值逐列导出**（未做重命名、未做 U/T 转换）。
- **规模**：644 条记录 = 基因组上的 tRNA 基因拷贝（同一 unique 序列可对应多个拷贝）。
  按 `Chr` 分：染色体 585 条、chloroplast 37 条、mitochondrion 22 条。
- **列**（24 列）：`Database-ID, Chr, Start, Stop, AA, AC` + 17 个结构区段列 + `Whole_sequence`。
- **区段列**（5′→3′，共 17 段）：`p-1, Accstem, 8-9, Dstem, Dloop, Dstem, p26, Acstem, Acloop, Acstem, Vregion, Tstem, Tloop, Tstem, Accstem, p73, CCA`
  - 注意 `Dstem / Acstem / Tstem` 各出现 **两次**（列名重复），需按列位置区分：前一个为 5′ 半、后一个为 3′ 半。
  - `p-1` 即额外 5′ 碱基位（下游表二/图里叫 `p1`）；全表仅 His 的 12 条基因组记录有值，恒为 `G`。
- **字母表**：区段列与 `Whole_sequence` 均为 RNA（U）。
- **自洽性**：17 段直接拼接 == `Whole_sequence`，644/644 实测一致。

## 表二：plantRNA_database_unique.tsv

- **构造**：表一按 `Whole_sequence` 去重（RNA/DNA 归一后按完全相同序列归并）→ **644 → 220**，恰好等于参考里的 220 条 unique tRNA。
- **Database-ID**：保留该 unique 序列在表一中**首次出现**的那条记录（220/220 对应）。
- **ID 重命名**：`tRNA-<分室>-<AA>-<反密码子>-<序号>`，如 `tRNA-Nu-Asp-GTC-1`；
  `Chr` 列由基因组位置改为分室 `Nu / Mt / Pt`。
- **规模**：220 条 = **Nu 172 + Mt 18 + Pt 30**（对应 46 / 17 / 30 个 anticodon family）。
- **区段列改名**：`p1, Accstem_5p, p8-9, Dstem_5p, Dloop, Dstem_3p, p26, Acstem_5p, Acloop, Acstem_3p, Vregion, Tstem_5p, Tloop, Tstem_3p, Accstem_3p, p73, CCA`
  （不再有重名列）。
- **字母表（重要）**：区段列是 **RNA（U）**，但 `Whole_sequence` 是 **DNA（T）**；二者 `U→T` 后完全一致，220/220 实测。
- **`fasta_record` 列**：形如 `>tRNA-Nu-Pro-TGG-1\nGGGCA...` 的单条 FASTA 字符串，其中 `\n` 是**字面两字符**，不是真换行。
- **与比对参考一致**：220/220 的 `ID` 与序列和
  `/home/HL/reference/Ath/TAIR/tRNA-unique_sn_snoRNA.fa` 中的 tRNA 记录**完全相同**（该文件 589 条 = 220 tRNA + 287 snoRNA + 82 snRNA）。

## 与下游图表的关系

- 表二 = **完整 220 条**的 Sprinzl 区段注释，可覆盖参考里的全部 tRNA 拷贝。
- 下载原件中另有一个 sheet「93 个 anticodon 家族代表序列」，是从 220 条里**每个 family 取一条**的产物
  （人工策展，挑选规则无文档记录）；其内容仍保留在
  `archive/20260902/paper/tRNA_dependent/dependent_site_overview/All-tRNA_rep-structure.csv`。
  下游图组已改用 `tRNA/count/ncRNA/family_top/family_top_tRNA.tsv` 的**表达量代表**（mean_CPM 最高），
  且位点按 family 聚合、220 条拷贝全部覆盖。
- 220 条中与所属 family 代表**长度不同**的只有 4 条，坐标投影需单独处理：
  `tRNA-Nu-Asn-GTT-5`（+1）、`tRNA-Nu-Leu-TAA-3`（+1）、`tRNA-Nu-Ser-AGA-5`（+2）、`tRNA-Mt-Tyr-GTA-1`（−8）；
  其余 216 条与代表等长且可无 gap 比对（纯碱基替换），区段边界可直接沿用。

## 脚本：add_sprinzl_pos.py

给任意位点表追加一列 `Sprinzl_pos`（canonical 0–76），放在 `pos` 列之后。

```bash
# 默认 chr / pos 列、默认注释（同目录表二）、默认输出 <stem>.Sprinzl<ext>
python3 add_sprinzl_pos.py -i tRNA/STAT/fisher/mRNA/trua1a_1_Flower_sites.positive.tsv
python3 add_sprinzl_pos.py -i tRNA/STAT/dependent_methylKit/ncRNA/Flower/trua1a/trua1a-Flower.dependent.tsv

# 列名 / 输出自定义
python3 add_sprinzl_pos.py -i in.tsv --chr-col chr --pos-col pos -o out.tsv
```

- **映射规则**与 `Figures/dependent_methylKit/tRNA/build_genotype_dependent_tables.sh` 完全一致：
  按区段累加定位 `region_pos`，再 `round_half_up((region_pos-1)*(L-1)/(n-1)) + 1` 到固定槽位；
  区段偏移合计 77 列（canonical 0–76）。
- **按拷贝注释**：直接用表二里该 `chr` 自己的区段序列，不借用 family 代表。
  220 条中有 8 条拷贝的区段长度与代表不同（D 环/茎的拆分不同或总长不同），用代表布局会把位点算偏；
  实测 468 条 dependent 记录中有 9 行落在这些拷贝上，其中 4 行的映射与代表布局不同
  （`Nu-Leu-TAA-3` pos 40 → 39、`Nu-Ser-AGA-5` pos 41 → 39），其余两者一致。
- **非 tRNA 行**（mRNA 转录本、sn/snoRNA 等表里没有的 `chr`）原样保留，`Sprinzl_pos` 留空；
  运行结束在 stderr 打印 `mapped / unknown / invalid pos` 统计与示例。
- `--strict`：遇到非法 `pos` 或越界时直接报错退出（默认只告警并留空）。
- 半值取上（`int(x+0.5)`），与 awk 版一致；不要换成 Python `round()`（银行家舍入）。

## 复核

```python
import csv

def read_tsv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]

_, t1 = read_tsv("plantRNA_database.tsv")          # 表一
_, t2 = read_tsv("plantRNA_database_unique.tsv")   # 表二
assert len(t1) == 644 and len(t2) == 220

# 表一：17 段拼接 == Whole_sequence（均为 RNA）
assert all("".join(row[6:23]) == row[23] for row in t1)

# 表二：17 段拼接 U->T == Whole_sequence（DNA）
assert all("".join(row[5:22]).replace("U", "T") == row[22] for row in t2)

# 表一去重序列集合 == 表二序列集合 == 220
assert {r[23].replace("U", "T") for r in t1} == {r[22] for r in t2}
```
