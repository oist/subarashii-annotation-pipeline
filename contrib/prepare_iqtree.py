#!/usr/bin/env python3
"""
Create IQTree/<gfX>/ folders for top-300 families by speciesrax_score
and copy their .trim files there.

Usage
    python prepare_iqtree.py  [inflation] [subset]

Defaults
    inflation = 1.8
    subset    = Normal
"""

import sys, subprocess, pathlib, shutil, pandas as pd

infl   = sys.argv[1] if len(sys.argv) > 1 else "1.8"
subset = sys.argv[2] if len(sys.argv) > 2 else "Normal"

# -------- locate result_dir via read_conf.sh ---------------------
SCRIPTS = pathlib.Path(__file__).resolve().parent
read_conf = SCRIPTS / "read_conf.sh"
conf_file = SCRIPTS / "pipeline.conf"

def conf_get(sec, key):
    cmd = ["bash", "-c", f"source {read_conf} && conf_get {conf_file} {sec} {key}"]
    return subprocess.check_output(cmd, text=True).strip()

RES = pathlib.Path(conf_get("dirs", "result_dir"))
NORMAL_DIR = RES / f"families/{infl}/{subset}"
SCORE_TSV  = RES / f"families/{infl}/family_scores.tsv"
IQ_DIR     = RES / f"families/{infl}/IQTree"
LOG_DIR    = RES / "logs"
IQ_DIR.mkdir(parents=True, exist_ok=True)
LOG_DIR.mkdir(exist_ok=True)

# -------- pick top 300 ------------------------------------------
df = pd.read_csv(SCORE_TSV, sep="\t")
top300 = df.sort_values("speciesrax_score", ascending=False).head(300)["family"]

filelist = LOG_DIR / f"iqtree_filelist_{infl}.txt"
with filelist.open("w") as lf:
    for fam in top300:
        src = NORMAL_DIR / f"{fam}.trim"
        if not src.exists():
            print(f"[WARN] missing {src}")
            continue
        dst_dir = IQ_DIR / fam
        dst_dir.mkdir(exist_ok=True)
        dst = dst_dir / src.name
        if not dst.exists():
            shutil.copy(src, dst)
        lf.write(str(dst) + "\n")

print(f"[prepare_iqtree] wrote {filelist} with {sum(1 for _ in open(filelist))} alignments")

