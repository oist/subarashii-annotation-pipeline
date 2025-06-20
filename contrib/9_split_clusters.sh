#!/usr/bin/env bash
# ------------------------------------------------------------
# 9_split_clusters.sh
# Wrapper that calls split_clusters.py with config-driven paths
# ------------------------------------------------------------
set -euo pipefail

CONF="pipeline.conf"
[ -f pipeline.local.conf ] && CONF="pipeline.local.conf"

source read_conf.sh

echo "[split_clusters] using config $CONF"
python split_clusters.py

