#!/bin/bash
# What-if: WT positive sites additionally require Mock gap ratio < 5% in every replicate
set -euo pipefail
cd "$(dirname "$0")"

# repo root = mock_gap5pct_check -> wildtype_landscape -> Figures -> PUSdependent
REPO_ROOT="$(cd ../../.. && pwd)"

/home/HL/miniforge3/envs/DeepPSU/bin/python mock_gap5pct_check.py \
    --repo-root "$REPO_ROOT" \
    --outdir .
