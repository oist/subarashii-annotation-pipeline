#!/usr/bin/env bash
# submit_step.sh  <step>     (e.g. eggnog)
set -euo pipefail

STEP="$1"

# -------- locate pipeline.conf wherever it lives ----------

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)   # absolute contrib/
ROOT_DIR=$(dirname "$SCRIPT_DIR")         # absolute repo root


# search order: CLI override > contrib/ > repo root
CONF="${PIPELINE_CONF:-}"
if [[ -z "$CONF" ]]; then
    [[ -f "$SCRIPT_DIR/pipeline.conf" ]] && CONF="$SCRIPT_DIR/pipeline.conf"
fi
[[ -z "$CONF" && -f "$ROOT_DIR/pipeline.conf" ]] && CONF="$ROOT_DIR/pipeline.conf"

[ -f "$CONF" ] || { echo "[submit] ERROR: cannot find pipeline.conf"; exit 1; }

[[ -f "$(dirname "$CONF")/pipeline.local.conf" ]] && LOCAL_CONF="$(dirname "$CONF")/pipeline.local.conf"
source "$SCRIPT_DIR/read_conf.sh"

ENGINE=$(conf_get "$CONF" cluster engine)
TEMPLATE_DIR=$(conf_get "$CONF" cluster template_dir)

if [[ "$ENGINE" == "none" ]]; then
    echo "[submit] Running $STEP locally ..."
    bash "$SCRIPT_DIR/${STEP}_run_local.sh"
    exit
fi

# -------- only SLURM for now --------
TMPL="$ROOT_DIR/$TEMPLATE_DIR/$ENGINE/${STEP}.sbatch.tmpl"
[ -f "$TMPL" ] || { echo "[submit] ERROR: template $TMPL not found"; exit 1; }

RES_DIR=$(conf_get "$CONF" dirs result_dir)
mkdir -p "$RES_DIR/logs"
SB="$RES_DIR/logs/${STEP}_$(date +%s).sbatch"

# pull resource values for this step
THREADS=$(conf_get "$CONF" $STEP threads)
MEM=$(conf_get "$CONF" $STEP mem)
TIME=$(conf_get "$CONF" $STEP walltime)

# substitute placeholders and write SBATCH file
# substitute placeholders and write SBATCH file

# substitute placeholders and write SBATCH file
sed -e "s|{threads}|$THREADS|g" \
    -e "s|{mem}|$MEM|g" \
    -e "s|{walltime}|$TIME|g" \
    -e "s|{root}|$ROOT_DIR|g" \
    -e "s|{conf}|$CONF|g" \
    -e "s|{logdir}|$RES_DIR/logs|g" \
    -e "s|{module}|$(conf_get "$CONF" $STEP module)|g" \
    -e "s|{filelist}|$MAFFT_FILELIST|g" \
    -e "s|{maxidx}|$MAFFT_MAXIDX|g" \
    "$TMPL" > "$SB"


echo "[submit] Submitting $SB"
sbatch "$SB"

