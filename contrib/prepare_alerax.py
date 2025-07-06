#!/usr/bin/env python3
"""
Build families_file.txt for AleRax from the top-300 IQ-TREE .ufboot files.

Usage
    python prepare_alerax.py [inflation]

Default
    inflation = 1.8
"""

import sys, subprocess, pathlib, shutil

infl = sys.argv[1] if len(sys.argv) > 1 else "1.8"

# ------------------------------------------------ locate dirs ----
SCRIPTS   = pathlib.Path(__file__).resolve().parent
read_conf = SCRIPTS / "read_conf.sh"
conf_file = SCRIPTS / "pipeline.conf"

def conf_get(sec, key):
    cmd = ["bash", "-c", f"source {read_conf} && conf_get {conf_file} {sec} {key}"]
    return subprocess.check_output(cmd, text=True).strip()

RES_DIR   = pathlib.Path(conf_get("dirs", "result_dir"))
IQ_DIR    = RES_DIR / f"families/{infl}/IQTree"
LOG_DIR   = RES_DIR / "logs"
LOG_DIR.mkdir(exist_ok=True)

families_file = RES_DIR / f"families/{infl}/IQTree/families_file.txt"

# ------------------------------------------------ collect ufboot --
ufboot_paths = list(IQ_DIR.glob("gf*/gf*.trim.ufboot"))
if len(ufboot_paths) != 300:
    print(f"[WARN] Expected 300 ufboot files, found {len(ufboot_paths)}")

with families_file.open("w") as out:
    out.write("[FAMILIES]\n")
    for p in sorted(ufboot_paths):
        fam = p.stem.split(".")[0]          # gf0123
        out.write(f"- {fam}\n")
        out.write(f"gene_tree = {p.resolve()}\n")

print(f"[prepare_alerax] wrote {families_file}")

