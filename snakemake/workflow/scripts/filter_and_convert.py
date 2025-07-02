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
import csv

e_cut   = float(snakemake.params["evalue_cut"])
cov_cut = float(snakemake.params["cov_cut"])

in_tsv  = snakemake.input[0]
out_abc = snakemake.output[0]

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

