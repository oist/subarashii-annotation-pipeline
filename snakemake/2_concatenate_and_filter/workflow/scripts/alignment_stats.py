#!/usr/bin/env python3
"""
Compute statistics for every *.trim file after BMGE/trimAl.

Output (tab-separated):
  family              gfX (from filename)
  nseqs               number of sequences
  longest             longest gap-stripped length
  shortest            shortest gap-stripped length
  median              median length
  paralog_in_species  Yes/No  (≥1 species has >1 gene)

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-09)
"""
import argparse
import glob
import sys, subprocess, pathlib, statistics
from collections import Counter
from Bio import SeqIO

# ── helpers ───────────────────────────────────────────────────────
def ungapped_len(record):
    return len(str(record.seq).replace("-", ""))

# ── iterate files ────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(descritpiton=__doc__, formatter_class=argparse.RawHelpTextFormatter)
    parser.add_argument("input", help="input directory with msa alignments or paths to msa alignemnts to be included", nargs="+")
    parser.add_argument("-o", "--output", help="output file to save statistics in")
    args = parser.parse_args()

    trim_files = args.input
    if len(args.input) == 1:
        trim_files = sorted(glob.glob("{}/*".format(args.input)))
    if not trim_files:
        sys.exit(f"[error] no .trim files in {align_dir}")

    rows = []
    for trim in trim_files:
        fam = trim.stem.split(".")[0]              # gf123 from gf123.trim
        recs = list(SeqIO.parse(trim, "fasta"))
        if not recs:
            continue

        lengths = [ungapped_len(r) for r in recs]

        # paralog test identical to split_clusters.py
        sp_counts = Counter(r.id.split("_g")[0] for r in recs)
        has_paralog = any(v > 1 for v in sp_counts.values())

        rows.append((
            fam,
            len(recs),
            max(lengths),
            min(lengths),
            int(statistics.median(lengths)),
            "Yes" if has_paralog else "No"
        ))

    with open(args.output, "w") as f:
        f.write("family\tnseqs\tlongest\tshortest\tmedian\tparalog_in_species\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")

    print(f"[align-stats] wrote {out_tsv} ({len(rows)} families)")


if __name__ == "__main__":
    main()
    exit(0)
