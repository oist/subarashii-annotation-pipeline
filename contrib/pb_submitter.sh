#!/usr/bin/env bash
# Build file list of *.phy and submit a SLURM array on Euler.
set -euo pipefail

SCRIPT_DIR=$(dirname "$0")
CONF="$SCRIPT_DIR/pipeline.conf"
source "$SCRIPT_DIR/read_conf.sh"

INFL=1.8                    # change if you run another inflation
PB_DIR=$(conf_get "$CONF" dirs result_dir)/families/${INFL}/PhyloBayes
LOG_DIR=$(conf_get "$CONF" dirs result_dir)/logs
mkdir -p "$LOG_DIR"

FILELIST="$LOG_DIR/pb_filelist_${INFL}.txt"
find "$PB_DIR" -name '*.phy' | sort > "$FILELIST"
TOTAL=$(wc -l < "$FILELIST")
(( TOTAL > 0 )) || { echo "[pb-submit] no .phy files found"; exit 1; }

# 1 task = 1 family
export PB_FILELIST="$FILELIST"
export PB_MAXIDX=$(( TOTAL - 1 ))

# submit through generic wrapper
bash "$SCRIPT_DIR/submit_step.sh" pb

