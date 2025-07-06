#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(dirname "$0")
CONF="$SCRIPT_DIR/pipeline.conf"
source "$SCRIPT_DIR/read_conf.sh"

INFL=1.8
#python "$SCRIPT_DIR/prepare_alerax.py" "$INFL"

RES=$(conf_get "$CONF" dirs result_dir)
FAMILIES_FILE="$RES/families/${INFL}/IQTree/families_file.txt"
OUT_PREFIX="$RES/families/${INFL}/IQTree/alerax_out"

export ALERAX_FFILE="$FAMILIES_FILE"
export ALERAX_PREFIX="$OUT_PREFIX"

bash "$SCRIPT_DIR/submit_step.sh" alerax

