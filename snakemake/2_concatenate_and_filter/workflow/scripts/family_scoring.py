#!/usr/bin/env python3
"""
Rank trimmed alignments for concatenation markers
and SpeciesRax gene-tree inference.

Output:
    results/families/<inflation>/family_scores.tsv  (tab-separated)
Columns:
    family, n_species, n_total, median_len, single_copy,
    dist_to_2, oversize_penalty,
    concat_score, speciesrax_score

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-14)
"""
import argparse
import glob
import sys, subprocess, pathlib, statistics, math
from collections import Counter, defaultdict
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawHelpTextFormatter)
    parser.add_argument("input", help="Direcotry containing msas or list of msa files to evaluate", nargs="+")
    parser.add_argument("-o", "--output", help="Output tsv file", default="family_scores.tsv")
    args = parser.parse_args()

    trim_path = args.input
    if len(args.input) == 1:
        trim_paths = sorted(glob.glob("{}/*".format(args.input)))
    if not trim_paths:
        sys.exit(f"[error] no .trim files in {ALIGN_DIR}")

    # -------- first pass: collect species universe & per-family data -
    total_species_set = set()
    fam_data = {}

    def ungapped(rec):
        return len(str(rec.seq).replace("-", ""))

    for trim in trim_paths:
        fam = trim.stem.split(".")[0]          # gf123 from gf123.trim
        recs = list(SeqIO.parse(trim, "fasta"))
        if not recs:
            continue

        counts = Counter(r.id.split("_g")[0] for r in recs)
        total_species_set.update(counts.keys())

        lengths = [ungapped(r) for r in recs]

        fam_data[fam] = dict(
            n_species=len(counts),
            n_total=len(recs),
            median_len=int(statistics.median(lengths)),
            single_copy=1 if all(v == 1 for v in counts.values()) else 0,
            dist_to_2=math.sqrt(sum((v - 2) ** 2 for v in counts.values())),
            counts=counts                          # keep for oversize calc
        )

    TOTAL_SPECIES = len(total_species_set)
    if TOTAL_SPECIES == 0:
        sys.exit("[error] could not detect any species codes")

    # -------- scoring ------------------------------------------------
    rows = []
    for fam, d in fam_data.items():
        oversize_penalty = max(0, d["n_total"] - 3 * d["n_species"]) / d["n_species"]
        universality_fraction = d["n_species"] / TOTAL_SPECIES

        concat_score = d["single_copy"] * (10 * d["n_species"] + 0.01 * d["median_len"])

        speciesrax_score = (
            universality_fraction * (0.01 * d["median_len"] + d["n_species"])
        ) / (1 + d["dist_to_2"] + oversize_penalty)

        rows.append((
            fam,
            d["n_species"], d["n_total"], d["median_len"],
            d["single_copy"],
            round(d["dist_to_2"], 2),
            round(oversize_penalty, 2),
            round(concat_score, 2),
            round(speciesrax_score, 4)
        ))

    # -------- write table -------------------------------------------
    with open(args.output, "w") as f:
        f.write("\t".join([
            "family", "n_species", "n_total", "median_len",
            "single_copy", "dist_to_2", "oversize_penalty",
            "concat_score", "speciesrax_score"]) + "\n")
        for r in sorted(rows, key=lambda x: x[-1], reverse=True):
            f.write("\t".join(map(str, r)) + "\n")

    print(f"[family_scoring] wrote {OUT_TSV} ({len(rows)} families)")

if __name__ == "__main__":
    main()
    exit(0)

