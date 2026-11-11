# BID-transcript Snakemake Configuration Guide

Two workflows (located in the `workflows/` directory) share the following preprocessing configuration:

- `workflows/polyA_mapping.smk`
- `workflows/totalRNA_mapping.smk`

## 1. Paired-end Data

```yaml
pipeline:
  read_layout: paired   # Can also be set to PE
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

Paired-end mode requires both `read1` and `read2`. The polyA workflow
preprocesses both reads and then uses R2 for all downstream alignment; the
totalRNA workflow preprocesses and analyzes R2 only.

## 2. Single-end Data

```yaml
pipeline:
  read_layout: single   # Can also be set to SE
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

For single-end FASTQ data, provide `read1` only; `read2` is not required.
Trimming, UMI extraction, alignment, and log statistics all use R1.

## 3. Deduplication Toggle

```yaml
pipeline:
  rmdup: true
```

- `true`: Perform deduplication. polyA uses `fastp --dedup`, totalRNA uses `seqkit rmdup`.
- `false`: Do not generate `*.dup.fq.gz`; trimmed FASTQ files pass directly to `umi_extract`.

## 4. UMI Region Trimming

`trim_front` and `trim_tail` must be non-negative integers giving the number
of bases trimmed from the front and tail of the read, respectively (as in
fastp).

To process only R1 in paired-end data, set R2 to 0:

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

To process only R2, reverse the settings. When the R1/R2 blocks or individual
trim fields are omitted, the per-workflow defaults apply:

| Workflow | R1 front/tail | R2 front/tail |
|---|---:|---:|
| polyA | 6 / 21 | 21 / 6 |
| totalRNA | 17 / 0 | 17 / 0 |

## 5. trim_galore Minimum Length

`trim_galore --length` controls the minimum read length (bp) retained after
trimming. Set it via `pipeline.trim_galore_length`, which must be a
non-negative integer. When omitted, the per-workflow default applies:

| Workflow | Default |
|---|---:|
| polyA | 82 |
| totalRNA | 35 |

## 6. Optional Expression Counting

`polyA_mapping.smk` and `totalRNA_mapping.smk` support a
`pipeline.expression` block. The feature is disabled when the block is
omitted, and the `rule all` targets are unchanged.

```yaml
pipeline:
  ref:
    rRNA: /path/to/ath_rRNA.fa
    mRNA: /path/to/cDNA.fa
  expression:
    enabled: true
    count_mode:
      rRNA: all_alignments   # Typically ncRNA for totalRNA projects
      mRNA: primary
```

- `count_mode` keys must match `pipeline.ref` exactly.
- `all_alignments`: Counts all mapped alignments in the filtered BAM; a read can count towards multiple conserved rRNA/ncRNA features.
- `primary`: Skips secondary/supplementary alignments; a read counts toward at most one mRNA.
- Two matrices are output for each `ref_type` (raw integer counts, no TPM/CPM normalization):
  - `count/<ref_type>/Mock.<mode>_counts.tsv`: Mock samples only, retaining the existing structure and naming.
  - `count/<ref_type>/all_samples.<mode>_counts.tsv`: All samples (mock and treat), in the order specified in config.yaml.

## 7. Compatibility with Existing Configurations

When a configuration omits the expression-count fields,
`polyA_mapping.smk` and `totalRNA_mapping.smk` disable the feature and fall
back to the default preprocessing configuration:

```yaml
pipeline:
  read_layout: paired
  rmdup: true
```

UMI trimming then uses the per-workflow defaults listed in the table above,
so existing paired-end configurations continue to work unchanged.
