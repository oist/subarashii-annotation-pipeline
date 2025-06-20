#!/usr/bin/env python3
"""
Read DIAMOND hits (tab-separated), keep only matches that pass the
e-value and coverage cut-offs defined in pipeline.conf, and write the
result in MCL .abc format (q s weight).

Inputs  (from config [dirs]/result_dir):
    all_vs_all.tsv     – produced by 6_diamond_all_vs_all.sh

Outputs (same directory):
    all_vs_all.abc
"""


import configparser
import pathlib
import csv

# ---------- read config (local overrides global) ----------
cfg = configparser.ConfigParser(inline_comment_prefixes=(";", "#"))   #  ← NEW
cfg.read(["pipeline.conf", "pipeline.local.conf"])                    #  ← unchanged

res_dir = pathlib.Path(cfg["dirs"]["result_dir"])
e_cut   = float(cfg["diamond"]["evalue_cut"])
cov_cut = float(cfg["diamond"]["cov_cut"])

in_tsv  = res_dir / "all_vs_all.tsv"
out_abc = res_dir / "all_vs_all.abc"

kept = 0
total = 0

with in_tsv.open() as fi, out_abc.open("w") as fo:
    reader = csv.reader(fi, delimiter="\t")
    for q, s, pid, bits, e, qcov, scov in reader:
        total += 1
        if float(e) > e_cut:
            continue
        if min(float(qcov), float(scov)) < cov_cut:
            continue
        fo.write(f"{q} {s} {bits}\n")
        kept += 1

print(f"[filter_and_convert] kept {kept}/{total} hits "
      f"(E <= {e_cut}, cov >= {cov_cut}%)")

