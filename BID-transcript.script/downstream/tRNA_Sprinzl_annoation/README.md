# tRNA Annotation and Mapping

Annotation databases and scripts for mapping native tRNA site coordinates to
canonical Sprinzl coordinates (0-76).

## Databases

- `plantRNA_database.tsv`: Original annotations from the plantRNA database (644 genomic copies).
- `plantRNA_database_unique.tsv`: Curated list of 220 unique tRNA sequences used as the primary reference for coordinate mapping.

## Scripts

### add_sprinzl_pos.py
Adds a `Sprinzl_pos` column (canonical 0-76) to a site table based on the reference databases.

**Usage:**

```bash
# Default mapping using native coordinate and database
python3 add_sprinzl_pos.py -i <input_file.tsv>

# Custom column mapping
python3 add_sprinzl_pos.py -i <input_file.tsv> --chr-col <chrom_col> --pos-col <pos_col> -o <output_file.tsv>
```

- **Mapping Rules**: Uses the tRNA-specific structural region annotations to normalize site positions into canonical 0-76 slots. Rows for sequences absent from the database are preserved with `Sprinzl_pos` left blank.
- **Strict Mode**: Use the `--strict` flag to exit with an error when a row has an invalid or out-of-range position.
