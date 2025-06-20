#!/usr/bin/env bash
set -euo pipefail
CONF="pipeline.conf"
source read_conf.sh

FA_DIR=$(conf_get "$CONF" dirs proteome_dir)
RES_DIR=$(conf_get "$CONF" dirs result_dir)

mkdir -p "$RES_DIR"
cat "$FA_DIR"/*.faa > "$RES_DIR/all_proteins.faa"
diamond makedb --in "$RES_DIR/all_proteins.faa" \
               -d  "$RES_DIR/proteomes.dmnd"

