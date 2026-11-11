#!/usr/bin/env bash
# Build WT tissue ratio distribution and tissue-pair Welch t-tests.

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

Rscript "${script_dir}/plot.R"
