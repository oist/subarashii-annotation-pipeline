#!/usr/bin/env python3
"""
Filter alignments by universality and copy the alignments passing the test into the provided output folder.

Universality = number of species in this gene family / maximum number of species in gene families

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-17)
"""
import argparse
import os, glob, pathlib, shutil
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", help="Directory with fasta files", required=True, type=pathlib.Path)
    parser.add_argument("-o", "--output", help="Output directory base to copy alignments to", required=True, type=pathlib.Path)
    parser.add_argument("-t", "--threshold", help="Universality threshold to use", default=0.8, type=float)
    parser.add_argument("-m", "--metadata", help="Metadata tsv or csv file containing a per-family data", required=True, type=pathlib.Path)
    args = parser.parse_args()

    # load scores
    df = pd.read_csv(args.metadata, sep="\t")
    df["universality"] = df["n_species"] / df["n_species"].max()

    # select families passing the threshold
    concat_fams = df[df["universality"] >= args.threshold]["family"].tolist()
    print(f"[concat] {len(concat_fams)} families ≥{args.threshold} % universality")

    # copy files passing the threshold
    args.output.mkdir(parents=True, exist_ok=True)
    missing_families = []
    for fam in concat_fams:
        trim = glob.glob(f"{args.input}/{fam}*")
        if len(trim) == 0:
            missing_families.append(fam)
            continue
        shutil.copy(trim[0], "{}/{}".format(args.output, os.path.basename(trim[0])))

    if len(missing_families) > 1:
        print("The following families' sequence file could not be found:")
        for f in missing_families:
            print(f)

if __name__ == "__main__":
    main()
    exit(0)
