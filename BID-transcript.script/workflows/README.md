# Workflows

Snakemake workflows for BID-seq library mapping and processing. Both workflows
analyze R2 only: trim, dedup/UMI clip, bowtie2 against the primary ncRNA or
rRNA reference, STAR against the secondary cDNA reference, then samFilter,
realignGap, and adjustGap.

## Available Workflows

- `polyA_mapping.smk`: Pipeline for poly(A)-enriched RNA libraries.
- `totalRNA_mapping.smk`: Pipeline for total RNA libraries.

## Usage

To run a workflow, use the `snakemake` command:

```bash
snakemake -s <workflow.smk> --configfile config.yaml --cores <num_cores>
```

The `config.yaml` file should be derived from `config.example.yaml`. Refer to `CONFIG_GUIDE.md` for details on configuration fields and settings.
