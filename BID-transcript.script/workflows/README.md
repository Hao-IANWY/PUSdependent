# workflows/ — Snakemake 工作流

| 文件 | 定位 |
|---|---|
| `polyA_mapping.smk` | polyA 主流程（含分块 adjustGap） |
| `totalRNA_mapping.smk` | totalRNA 流程 |

## 运行

```bash
snakemake -s workflows/polyA_mapping.smk --configfile config.yaml --cores 120
```

config.yaml 从 [`../config.example.yaml`](../config.example.yaml) 复制修改，
字段说明（单/双端、UMI、去重、trim_galore 长度等）见 [`../CONFIG_GUIDE.md`](../CONFIG_GUIDE.md)。

smk 不硬编码任何工具的相对位置：流程内工具通过 config 的 `paths.script_dir`
引用，移动本仓库位置只需改 config 里的路径。

## 表达 count 输出

`polyA_mapping.smk` 与 `totalRNA_mapping.smk` 支持可选 `pipeline.expression` 配置。启用后每个 ref_type 固定生成两张原始 count 矩阵：

- `all_alignments`：统计所有 mapped alignment，适用于保守 rRNA/ncRNA，一条 read 可进入多个 feature；
- `primary`：只统计 primary alignment，适用于 mRNA，一条 read 至多进入一个 transcript；
- `count/<ref_type>/Mock.<mode>_counts.tsv`：仅 mock 组样品，沿用既有结构；
- `count/<ref_type>/all_samples.<mode>_counts.tsv`：全部样品（mock 与 treat）。

矩阵包含参考 FASTA 中的全部 feature 和原始整数 count。

两类模式的计数总和语义不同，不能直接相加或直接比较总量。
