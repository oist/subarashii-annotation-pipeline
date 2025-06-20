#!/usr/bin/env bash
# ------------------------------------------------------------
# 10_run_eggnog.sh
# Run EggNOG-mapper on all renamed proteins in one file.
# ------------------------------------------------------------
set -euo pipefail

CONF="pipeline.conf"
[ -f pipeline.local.conf ] && CONF="pipeline.local.conf"

source read_conf.sh        # conf_get()

# ----- pull config -----
RES_DIR=$(conf_get "$CONF" dirs result_dir)
FASTA="$RES_DIR/all_proteins.faa"
OUT_DIR="$RES_DIR/eggnog"
DB_DIR=$(conf_get "$CONF" eggnog db_dir)
THREADS=$(conf_get "$CONF" eggnog threads)
MODE=$(conf_get "$CONF" eggnog mode)

mkdir -p "$OUT_DIR"

echo "[eggnog] running emapper on $(basename "$FASTA")"
emapper.py \
    -i "$FASTA" \
    --cpu "$THREADS" \
    --data_dir "$DB_DIR" \
    -m "$MODE" \
    --output "$OUT_DIR/annotations" \
    --override

echo "[eggnog] output files in $OUT_DIR/"

