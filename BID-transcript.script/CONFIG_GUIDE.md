# BID-transcript Snakemake 配置说明

三个工作流（位于 `workflows/` 目录）共用以下预处理配置：

- `workflows/polyA_mapping.smk`
- `workflows/totalRNA_mapping.smk`
- `workflows/polyA_mapping-adjustDeletion.smk`

## 1. 双端数据

```yaml
pipeline:
  read_layout: paired   # 也可以写 PE
  rmdup: true
  umi_extract:
    R1:
      trim_front: 6
      trim_tail: 21
    R2:
      trim_front: 21
      trim_tail: 6

samples:
  sampleA:
    sample_id: sampleA
    mock:
      rep1:
        file_index: Mock.sampleA.rep1
        read1: raw_data/Mock.sampleA.rep1.R1.fq.gz
        read2: raw_data/Mock.sampleA.rep1.R2.fq.gz
    treat:
      rep1:
        file_index: Treat.sampleA.rep1
        read1: raw_data/Treat.sampleA.rep1.R1.fq.gz
        read2: raw_data/Treat.sampleA.rep1.R2.fq.gz
```

双端模式必须同时提供 `read1` 和 `read2`。polyA 工作流会对两端进行预处理，并沿用原流程以 R2 进行后续比对；totalRNA 工作流沿用原流程，只分析 R2。

## 2. 单端数据

```yaml
pipeline:
  read_layout: single   # 也可以写 SE
  rmdup: false
  umi_extract:
    R1:
      trim_front: 8
      trim_tail: 0

samples:
  sampleA:
    sample_id: sampleA
    mock:
      rep1:
        file_index: Mock.sampleA.rep1
        read1: raw_data/Mock.sampleA.rep1.fq.gz
    treat:
      rep1:
        file_index: Treat.sampleA.rep1
        read1: raw_data/Treat.sampleA.rep1.fq.gz
```

单端 FASTQ 统一填写在 `read1`，不要求 `read2`。后续 trim、UMI 裁剪、比对和日志统计均使用 R1。

## 3. 去重开关

```yaml
pipeline:
  rmdup: true
```

- `true`：执行去重。polyA 使用 `fastp --dedup`，totalRNA 使用 `seqkit rmdup`。
- `false`：不生成 `*.dup.fq.gz`，trim 后的 FASTQ 直接进入 `umi_extract`。

## 4. UMI 区域裁剪

`trim_front` 和 `trim_tail` 必须是大于或等于 0 的整数，分别对应 fastp 的前端和末端裁剪碱基数。

如果双端数据只需要处理 R1，可以把 R2 设置为 0：

```yaml
pipeline:
  umi_extract:
    R1:
      trim_front: 8
      trim_tail: 4
    R2:
      trim_front: 0
      trim_tail: 0
```

只需要处理 R2 时，反向设置即可。未填写的 R1/R2 或单个裁剪字段会采用对应工作流的历史默认值：

| 工作流 | R1 front/tail | R2 front/tail |
|---|---:|---:|
| polyA / adjustDeletion | 6 / 21 | 21 / 6 |
| totalRNA | 17 / 0 | 17 / 0 |

## 5. trim_galore 最小长度

`trim_galore --length` 控制 trim 后保留的最小 read 长度（bp）。可以通过 `pipeline.trim_galore_length` 指定，必须是大于或等于 0 的整数。不填时使用对应工作流的历史默认值：

| 工作流 | 默认值 |
|---|---:|
| polyA / adjustDeletion | 82 |
| totalRNA | 35 |

## 6. 可选表达 count

`polyA_mapping.smk` 和 `totalRNA_mapping.smk` 支持 `pipeline.expression`。不填写该配置时功能关闭，旧 `rule all` 目标不变。

```yaml
pipeline:
  ref:
    rRNA: /path/to/ath_rRNA.fa
    mRNA: /path/to/cDNA.fa
  expression:
    enabled: true
    count_mode:
      rRNA: all_alignments   # totalRNA 项目此处通常为 ncRNA
      mRNA: primary
```

- `count_mode` 的 key 必须与 `pipeline.ref` 完全一致；
- `all_alignments` 统计 filtered BAM 中所有 mapped alignment，一条 read 可同时计入多个保守 rRNA/ncRNA feature；
- `primary` 跳过 secondary/supplementary alignment，一条 read 至多计入一个 mRNA；
- 每个 ref_type 固定输出两张矩阵，均为原始整数 count，不做 TPM/CPM：
  - `count/<ref_type>/Mock.<mode>_counts.tsv`：仅 mock 组样品，沿用既有结构与命名；
  - `count/<ref_type>/all_samples.<mode>_counts.tsv`：全部样品（mock 与 treat），样品顺序同 config.yaml。

## 7. 向后兼容

旧 YAML 不包含表达 count 字段时，`polyA_mapping.smk` 与 `totalRNA_mapping.smk` 默认关闭该功能；预处理配置默认：

```yaml
pipeline:
  read_layout: paired
  rmdup: true
```

UMI 裁剪使用上表中的历史默认值，因此旧双端配置可继续使用。
