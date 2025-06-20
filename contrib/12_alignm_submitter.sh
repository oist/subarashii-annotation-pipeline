#!/usr/bin/env bash
# Generates a file list & submits via submit_step.sh

module load stack/.2024-05-silent  gcc/13.2.0
module load mafft/7.505

set -euo pipefail
SCRIPT_DIR=$(dirname "$0")
CONF="$SCRIPT_DIR/pipeline.conf"          # adjust if needed
source "$SCRIPT_DIR/read_conf.sh"

RES_DIR=$(conf_get "$CONF" dirs result_dir)
INFL=$(conf_get "$CONF" mafft inflation)
NORM_DIR="$RES_DIR/families/${INFL}/Normal"
LOG_DIR=$(conf_get "$CONF" dirs log_dir)
mkdir -p "$LOG_DIR"

FILELIST="$LOG_DIR/mafft_filelist_${INFL}.txt"
ls "$NORM_DIR"/*.faa > "$FILELIST"

# tell the template how many tasks
MAXIDX=$(( $(wc -l < "$FILELIST") - 1 ))

export MAFFT_FILELIST="$FILELIST"
export MAFFT_MAXIDX="$MAXIDX"

bash "$SCRIPT_DIR/submit_step.sh" mafft

