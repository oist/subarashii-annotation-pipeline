#!/usr/bin/env bash
# ------------------------------------------------------------
# 7_filter_and_convert.sh
# Wrapper that calls the Python filter script so the pipeline
# numbering stays consistent.
# ------------------------------------------------------------
set -euo pipefail

# allow users to install custom virtual env, otherwise rely on PATH
python filter_and_convert.py

