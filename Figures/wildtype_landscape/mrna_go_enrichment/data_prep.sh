#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "${script_dir}/../../.." && pwd)"
input_dir="${project_root}/polyARNA/STAT/fisher/mRNA"

python3 - "$input_dir" "$script_dir" <<'PY'
from pathlib import Path
import csv
import sys

input_dir = Path(sys.argv[1])
out_dir = Path(sys.argv[2])

samples = {
    "Root": input_dir / "WT-Root_sites.positive.tsv",
    "Leaf": input_dir / "WT-leaf_sites.positive.tsv",
    "Flower": input_dir / "WT-Flower_sites.positive.tsv",
}

missing = [str(path) for path in samples.values() if not path.is_file()]
if missing:
    raise SystemExit("Missing WT mRNA positive-site files:\n" + "\n".join(missing))

site_rows = {}
site_tissues = {}
tissue_sites = {tissue: set() for tissue in samples}

for tissue, path in samples.items():
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"chr", "pos", "strand"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(f"Unexpected columns in {path}: {reader.fieldnames}")
        for row in reader:
            key = (row["chr"], int(row["pos"]), row["strand"])
            if tissue in site_tissues.get(key, set()):
                raise SystemExit(f"Duplicate site in {path}: {key}")
            site_rows.setdefault(key, row)
            site_tissues.setdefault(key, set()).add(tissue)
            tissue_sites[tissue].add(key)

root_keys = tissue_sites["Root"]
leaf_keys = tissue_sites["Leaf"]
flower_keys = tissue_sites["Flower"]

groups = {
    "root_unique": root_keys - leaf_keys - flower_keys,
    "leaf_unique": leaf_keys - root_keys - flower_keys,
    "flower_unique": flower_keys - root_keys - leaf_keys,
    "common": root_keys & leaf_keys & flower_keys,
}

# Pairwise-only categories are retained in the QC summary, but are not GO
# categories because the reference analysis contrasts all-three common sites
# with sites detected in only one tissue.
pairwise_counts = {
    "Root_Leaf_shared_not_Flower": len((root_keys & leaf_keys) - flower_keys),
    "Root_Flower_shared_not_Leaf": len((root_keys & flower_keys) - leaf_keys),
    "Leaf_Flower_shared_not_Root": len((leaf_keys & flower_keys) - root_keys),
}

gene_files = {
    "root_unique": out_dir / "root_unique_gene_list.tsv",
    "leaf_unique": out_dir / "leaf_unique_gene_list.tsv",
    "flower_unique": out_dir / "flower_unique_gene_list.tsv",
    "common": out_dir / "common_gene_list.tsv",
}

summary_rows = []
for tissue in samples:
    keys = tissue_sites[tissue]
    summary_rows.append({
        "category": f"{tissue}_input",
        "site_count": len(keys),
        "gene_count": len({key[0].split(".", 1)[0] for key in keys}),
    })

for group, keys in groups.items():
    genes = sorted({key[0].split(".", 1)[0] for key in keys})
    gene_files[group].write_text("\n".join(genes) + ("\n" if genes else ""))

    site_file = out_dir / f"{group}_sites.tsv"
    with site_file.open("w", newline="") as handle:
        fieldnames = ["chr", "pos", "strand", "ratio", "depth", "present_in"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for key in sorted(keys):
            row = dict(site_rows[key])
            row["pos"] = key[1]
            row["present_in"] = ";".join(
                tissue for tissue in samples if tissue in site_tissues[key]
            )
            writer.writerow(row)

    summary_rows.append({
        "category": group,
        "site_count": len(keys),
        "gene_count": len(genes),
    })

# All modified genes per tissue (every gene carrying >= 1 positive site in
# that tissue, regardless of sharing). These three lists feed the combined GO
# figure, whose groups are Roots / Leaves / Flowers.
tissue_gene_files = {
    "Root": out_dir / "root_all_gene_list.tsv",
    "Leaf": out_dir / "leaf_all_gene_list.tsv",
    "Flower": out_dir / "flower_all_gene_list.tsv",
}
for tissue, gene_file in tissue_gene_files.items():
    genes = sorted({key[0].split(".", 1)[0] for key in tissue_sites[tissue]})
    gene_file.write_text("\n".join(genes) + ("\n" if genes else ""))

for category, count in pairwise_counts.items():
    summary_rows.append({"category": category, "site_count": count, "gene_count": "NA"})

summary_file = out_dir / "mrna_go_input_summary.tsv"
with summary_file.open("w", newline="") as handle:
    fieldnames = ["category", "site_count", "gene_count"]
    writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(summary_rows)

print(f"[data_prep] Wrote gene lists and site tables to: {out_dir}")
print(
    "[data_prep] Input sites: "
    + ", ".join(f"{tissue}={len(tissue_sites[tissue])}" for tissue in samples)
)
print(
    "[data_prep] GO site categories: "
    + ", ".join(f"{group}={len(keys)}" for group, keys in groups.items())
)
print("[data_prep] Pairwise-only QC: " + ", ".join(f"{k}={v}" for k, v in pairwise_counts.items()))
print(
    "[data_prep] All-tissue gene lists: "
    + ", ".join(
        f"{tissue}={len({key[0].split('.', 1)[0] for key in tissue_sites[tissue]})}"
        for tissue in samples
    )
)
PY
