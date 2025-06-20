#!/usr/bin/env bash
# ------------------------------------------------------------
# 6_diamond_all_vs_all.sh
# Run an all-against-all DIAMOND search on the combined FASTA,
# writing tab-separated hits to results/all_vs_all.tsv.
# ------------------------------------------------------------
set -euo pipefail

# ---------- locate the active config file ----------
CONF="pipeline.conf"
[ -f pipeline.local.conf ] && CONF="pipeline.local.conf"

# ---------- tiny INI reader ----------
source read_conf.sh          # defines conf_get()

# ---------- pull settings ----------
RES_DIR=$(conf_get "$CONF" dirs result_dir)
THREADS=$(conf_get "$CONF" diamond threads)
DB="$RES_DIR/proteomes.dmnd"
QUERY="$RES_DIR/all_proteins.faa"
OUT="$RES_DIR/all_vs_all.tsv"

mkdir -p "$RES_DIR"

# ---------- run DIAMOND ----------
diamond blastp \
  --db       "$DB" \
  --query    "$QUERY" \
  --out      "$OUT" \
  --outfmt   6 qseqid sseqid pident bitscore evalue qcovhsp scovhsp \
  --more-sensitive \
  --threads  "$THREADS"

# ---------- guarantee self-hits so singletons survive ----------
# q s pident bits evalue qcov scov  (7 cols, tab-separated)
awk '/^>/{
        sub(/^>/,"",$0);
        printf "%s\t%s\t100\t1000\t1e-200\t100\t100\n",$0,$0
     }' "$QUERY" >> "$OUT"

