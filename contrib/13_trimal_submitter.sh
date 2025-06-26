#!/usr/bin/env bash
set -euo pipefail
set -x                                   # turn on tracing AFTER strict mode

SCRIPT_DIR=$(dirname "$0")
CONF="$SCRIPT_DIR/pipeline.conf"         # adjust if the real file is elsewhere
source "$SCRIPT_DIR/read_conf.sh"

echo "[debug] using conf: $CONF"

RESULTS=$(conf_get "$CONF" dirs result_dir)
INFL=$(conf_get "$CONF" families inflation)
ALIGN_DIR="$RESULTS/families/${INFL}/Normal"
echo "[debug] ALIGN_DIR=$ALIGN_DIR"

LOG_DIR="$RESULTS/logs"
mkdir -p "$LOG_DIR"

FILELIST="$LOG_DIR/trimal_filelist_${INFL}.txt"
echo "[debug] building file list → $FILELIST"
find "$ALIGN_DIR" -name '*.aln' | sort > "$FILELIST"
TOTAL=$(wc -l < "$FILELIST")
echo "[debug] TOTAL alignments = $TOTAL"

(( TOTAL > 0 )) || { echo "[error] No .aln files"; exit 1; }

MAX_JOBS=30
CHUNK=$(( (TOTAL + MAX_JOBS - 1) / MAX_JOBS ))
MAXIDX=$(( (TOTAL + CHUNK - 1) / CHUNK - 1 ))
export TRIMAL_FILELIST="$FILELIST" TRIMAL_CHUNK="$CHUNK" TRIMAL_MAXIDX="$MAXIDX"
echo "[debug] chunk=$CHUNK  array size=$((MAXIDX+1))"

bash "$SCRIPT_DIR/submit_step.sh" trimal

