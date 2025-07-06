#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(dirname "$0")
CONF="$SCRIPT_DIR/pipeline.conf"
source "$SCRIPT_DIR/read_conf.sh"

INFL=1.8                      # same inflation used for prepare_iqtree
#python "$SCRIPT_DIR/prepare_iqtree.py" "$INFL" Normal

FILELIST=$(conf_get "$CONF" dirs result_dir)/logs/iqtree_filelist_${INFL}.txt
MAXIDX=$(( $(wc -l < "$FILELIST") - 1 ))
export IQTREE_FILELIST="$FILELIST" IQTREE_MAXIDX="$MAXIDX"

bash "$SCRIPT_DIR/submit_step.sh" iqtree

