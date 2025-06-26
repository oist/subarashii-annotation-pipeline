#!/usr/bin/env bash
# Run on the login node or inside an interactive job.
set -euo pipefail
INFLATION=1.8
FSET=Normal
python alignment_stats.py "$INFLATION" "$FSET"

