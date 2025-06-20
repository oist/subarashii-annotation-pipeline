#!/usr/bin/env bash
# ------------------------------------------------------------
# 8_mcl_sweep.sh
# Cluster the all-vs-all similarity graph with MCL for every
# inflation value listed in pipeline.conf -> [mcl] inflation_vals
# ------------------------------------------------------------
set -euo pipefail

# ---- locate the active config file ----
CONF="pipeline.conf"
[ -f pipeline.local.conf ] && CONF="pipeline.local.conf"

# ---- tiny INI reader (defines conf_get) ----
source read_conf.sh

# ---- pull paths and parameters from the config ----
RES_DIR=$(conf_get "$CONF" dirs result_dir)          # e.g. results
ABC="$RES_DIR/all_vs_all.abc"                       # input for MCL
OUT_DIR="$RES_DIR/mcl"                              # output folder
INFLATIONS=$(conf_get "$CONF" mcl inflation_vals)   # "1.4 1.8 2.0"

mkdir -p "$OUT_DIR"

echo "[mcl_sweep] using $ABC"
echo "[mcl_sweep] inflation values: $INFLATIONS"

# ---- run MCL for every inflation value ----
for I in $INFLATIONS; do
    OUT_FILE="$OUT_DIR/clusters_I${I}.txt"
    echo "[mcl_sweep] running mcl -I $I  ->  $OUT_FILE"
    mcl "$ABC" --abc -I "$I" -o "$OUT_FILE"
done

echo "[mcl_sweep] done. Cluster files are in $OUT_DIR/"

