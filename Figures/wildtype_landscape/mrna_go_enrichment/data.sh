#!/usr/bin/env bash
# Entry point: build the GO-enrichment input tables via data_prep.sh.
set -euo pipefail
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "${script_dir}/data_prep.sh"
