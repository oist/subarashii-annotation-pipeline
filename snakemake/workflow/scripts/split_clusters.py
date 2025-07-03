#!/usr/bin/env python3
"""
Split an MCL cluster file into separate family FASTA files.

Creates subfolders in <result_dir>/families/<inflation>/:
   Singletons / Duplets / Triplets / Normal / Huge

Also writes two logs in the same directory:
   family_type_counts.tsv
   family_details.tsv

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-03)
"""

import argparse
import pathlib
import statistics
from collections import Counter, defaultdict
from Bio import SeqIO

# ------------------------- helpers -----------------------------
def family_folder(size):
    if size == 1:
        return "Singletons"
    if size == 2:
        return "Duplets"
    if size == 3:
        return "Triplets"
    if size >= huge_cutoff:
        return "Huge"
    return "Normal"

def lengths(records):
    return [len(r.seq) for r in records]

def has_paralog(records):
    species_seen = Counter(r.id.split("_g")[0] for r in records)
    return any(cnt > 1 for cnt in species_seen.values())

# ------------------------- main --------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", help="Output base directory. Subdirectory with the inflation parameter will be created inside of this directory", required=True)
    parser.add_argument("-i", "--inflation", help="Inlfation parameter", required=True)
    parser.add_argument("-f", "--fasta", help="Fasta formatted file containing all the sequences", required=True)
    parser.add_argument("-m", "--mcl", help="MCL's cluster output file", required=True)
    parser.add_argument("-c", "--cutoff", help="Cutoff number for huge families", required=True)
    args = parser.parse_args()

    res_dir = Path(args.output)
    infl = args.inflation
    huge_cutoff = args.cutoff

    cluster_file = args.mcl
    out_root = res_dir / infl
    out_root.mkdir(parents=True, exist_ok=True)

    seq_dict = {rec.id: rec for rec in SeqIO.parse(args.fasta, "fasta")}
    print(f"[split_clusters] loaded {len(seq_dict)} sequences")

    # iterate cluster
    type_counter = Counter()
    detail_rows  = []

    with args.mcl.open() as fh:
        for idx, line in enumerate(fh):
            ids = line.strip().split()
            size = len(ids)
            folder = family_folder(size)
            type_counter[folder] += 1

            gf_name = f"gf{idx}"
            out_dir = res_dir / folder
            out_dir.mkdir(exist_ok=True)
            out_faa = out_dir / f"{gf_name}.faa"

            recs = [seq_dict[i] for i in ids if i in seq_dict]
            SeqIO.write(recs, out_faa.open("w"), "fasta")

            lens = lengths(recs)
            detail_rows.append((
                gf_name, size,
                max(lens), min(lens),
                int(statistics.median(lens)),
                "Yes" if has_paralog(recs) else "No"
            ))

    # write logs
    with (out_root / "family_type_counts.tsv").open("w") as f:
        for t in ("Singletons", "Duplets", "Triplets", "Normal", "Huge"):
            f.write(f"{t}\t{type_counter[t]}\n")

    with (out_root / "family_details.tsv").open("w") as f:
        f.write("family\tsize\tlongest\tshortest\tmedian\tparalog_in_species\n")
        for row in detail_rows:
            f.write("\t".join(map(str, row)) + "\n")

    print("[split_clusters] done. Outputs in", out_root)

if __name__ == "__main__":
    main()
    exit(0)

