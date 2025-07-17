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

from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

def write_phylip_relaxed(records, outfile):
    aln = MultipleSeqAlignment(list(records))
    AlignIO.write(aln, outfile, "phylip-relaxed")

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", help="Directory with fasta files", required=True, type=pathlib.Path)
    parser.add_argument("-o", "--output", help="Output directory base to copy alignments to", required=True, type=pathlib.Path)
    parser.add_argument("-m", "--metadata", help="Metadata tsv or csv file containing a per-family data", required=True, type=pathlib.Path)
    args = parser.parse_args()

    # load scores
    df = pd.read_csv(args.metadata, sep="\t")

    # select top 300 families by speciesrax score
    top_pb = df.sort_values("speciesrax_score", ascending=False).head(300)

    args.output.mkdir(parents=True, exist_ok=True)
    missing_families = []
    for fam in top_pb["family"]:
        trim = glob.glob(f"{args.input}/{fam}*")
        if len(trim) == 0:
            missing_families.append(fam)
            continue
        recs = list(SeqIO.parse(trim[0], "fasta"))
        write_phylip_relaxed(recs, "{}/{}.phylip".format(args.output, ".".join(os.path.basename(trim[0]).split(".")[:-1])))

    print(f"[pb] prepared {len(top_pb)} families under {args.output}")

    if len(missing_families) > 1:
        print("The following families' sequence file could not be found:")
        for f in missing_families:
            print(f)

if __name__ == "__main__":
    main()
    exit(0)
