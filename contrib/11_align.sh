#!/usr/bin/env bash
# Align every .faa in the Normal/ folder locally (no SLURM).



set -euo pipefail
source read_conf.sh
CONF=${1:-contrib/pipeline.conf}

RES_DIR=$(conf_get "$CONF" dirs result_dir)
INFL=$(conf_get "$CONF" mafft inflation)
NORM_DIR="$RES_DIR/families/${INFL}/Normal"
THR=$(conf_get "$CONF" mafft threads)

echo "[mafft local] aligning in $NORM_DIR"
for faa in "$NORM_DIR"/*.faa; do
    [[ -f $faa ]] || continue
    mafft --auto --thread "$THR" "$faa" > "${faa%.faa}.aln"
done

