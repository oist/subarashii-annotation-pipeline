#!/usr/bin/env python3
"""
prepare_pb.py  –  select top families and set up folders for PhyloBayes.

Usage
  python 16_prepare_pb.py  [inflation] [family_set] [N]

Defaults
  inflation   = 1.8
  family_set  = Normal
  N           = 1000   (top-N families)

Creates:
  results/families/<infl>/PhyloBayes/gfX/gfX.trim  (symlink)
  plus a master file list: results/logs/pb_filelist_<infl>.txt
"""
import sys, subprocess, pathlib, shutil, pandas as pd

infl  = sys.argv[1] if len(sys.argv) > 1 else "1.8"
subset= sys.argv[2] if len(sys.argv) > 2 else "Normal"
TOPN  = int(sys.argv[3]) if len(sys.argv) > 3 else 1000

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
read_conf  = SCRIPT_DIR / "read_conf.sh"
conf_file  = SCRIPT_DIR / "pipeline.conf"

def conf_get(sec, key):
    cmd=["bash","-c",f"source {read_conf} && conf_get {conf_file} {sec} {key}"]
    return subprocess.check_output(cmd,text=True).strip()

RESULTS = pathlib.Path(conf_get("dirs","result_dir"))
ALIGN_DIR=RESULTS/f"families/{infl}/{subset}"
SCORE_TSV=RESULTS/f"families/{infl}/family_scores.tsv"
PB_ROOT  = RESULTS/f"families/{infl}/PhyloBayes"
LOG_DIR  = RESULTS/"logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)

df = pd.read_csv(SCORE_TSV, sep="\t")
df = df.sort_values("speciesrax_score", ascending=False).head(TOPN)

selected=[]
for fam in df["family"]:
    src = ALIGN_DIR/f"{fam}.trim"
    if not src.exists():
        print(f"[warn] missing {src}, skip"); continue
    dest_dir = PB_ROOT/fam
    dest_dir.mkdir(parents=True, exist_ok=True)
    link = dest_dir/f"{fam}.trim"
    if not link.exists():
        link.symlink_to(src)
    selected.append(link)

filelist = LOG_DIR/f"pb_filelist_{infl}.txt"
with filelist.open("w") as f:
    for p in selected: f.write(str(p)+"\n")

print(f"[prepare_pb] prepared {len(selected)} families, list → {filelist}")

