#!/usr/bin/env bash
# submit_step.sh  <step>     (e.g. eggnog)
set -euo pipefail

STEP="$1"

# -------- locate pipeline.conf wherever it lives ----------

SCRIPT_DIR=$(cd "$(dirname "$0")"; pwd)   # absolute contrib/
ROOT_DIR=$(dirname "$SCRIPT_DIR")         # absolute repo root

: "${IQTREE_FILELIST:=}" "${IQTREE_MAXIDX:=0}"
: "${ALERAX_FFILE:=}" "${ALERAX_PREFIX:=}"


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

# TMPL="$ROOT_DIR/$TEMPLATE_DIR/slurm/pb.sbatch.tmpl"



# ------------ resources (safe defaults if key missing) -------------
THREADS=$(conf_get "$CONF" "$STEP" threads 2>/dev/null || echo 1)
MEM=$(conf_get "$CONF" "$STEP" mem 2>/dev/null || echo 2G)
TIME=$(conf_get "$CONF" "$STEP" walltime 2>/dev/null || echo 01:00:00)

# ------------ optional placeholders ----------------
: "${TRIMAL_FILELIST:=}" "${TRIMAL_CHUNK:=}" "${TRIMAL_MAXIDX:=}"
: "${BMGE_FILELIST:=}"   "${BMGE_CHUNK:=}"   "${BMGE_MAXIDX:=}"
: "${JAR_PATH:=$(conf_get "$CONF" bmge jar_path 2>/dev/null || true)}"
: "${JAVA_MEM:=$(conf_get "$CONF" bmge java_mem 2>/dev/null || true)}"

SED_ARGS=(
  -e "s|{threads}|$THREADS|g"
  -e "s|{mem}|$MEM|g"
  -e "s|{walltime}|$TIME|g"
  -e "s|{root}|$ROOT_DIR|g"
  -e "s|{conf}|$CONF|g"
  -e "s|{logdir}|$RES_DIR/logs|g"
  -e "s|{module}|$(conf_get "$CONF" "$STEP" module 2>/dev/null || echo)|g"
)

case "$STEP" in
  trimal)
    SED_ARGS+=(
      -e "s|{filelist}|$TRIMAL_FILELIST|g"
      -e "s|{chunk}|$TRIMAL_CHUNK|g"
      -e "s|{maxidx}|$TRIMAL_MAXIDX|g"
    ) ;;
  bmge)
    SED_ARGS+=(
      -e "s|{filelist}|$BMGE_FILELIST|g"
      -e "s|{chunk}|$BMGE_CHUNK|g"
      -e "s|{maxidx}|$BMGE_MAXIDX|g"
      -e "s|{jar_path}|$JAR_PATH|g"
      -e "s|{java_mem}|$JAVA_MEM|g"
    ) ;;
    pb)
  SED_ARGS+=(
    -e "s|{filelist}|$PB_FILELIST|g"
    -e "s|{maxidx}|$PB_MAXIDX|g"
  ) ;;
  iqtree)
    SED_ARGS+=(
      -e "s|{filelist}|$IQTREE_FILELIST|g"
      -e "s|{maxidx}|$IQTREE_MAXIDX|g"
  ) ;;
  alerax)
    SED_ARGS+=(
      -e "s|{families_file}|$ALERAX_FFILE|g"
      -e "s|{out_prefix}|$ALERAX_PREFIX|g"
    ) ;;

esac

sed "${SED_ARGS[@]}" "$TMPL" > "$SB"
echo "[submit] wrote $SB"
sbatch "$SB"
echo "[submit] sbatch exit status $?"

