#!/usr/bin/env python3
"""
Prepare two alignment sets:

1. Concatenation candidates
   - universality >= 0.80
   - keep .trim
   - concatenate and write concatenated.faa + concatenated.phy

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-16)
"""
import argparse
import sys, glob, os, pathlib
from Bio import SeqIO

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def write_phylip_relaxed(records, outfile):
    """
    records : iterable of SeqRecord
    outfile : Path | str
    Writes PHYLIP-relaxed (seq type inferred).
    """
    aln = MultipleSeqAlignment(list(records))
    AlignIO.write(aln, outfile, "phylip-relaxed")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", help="Directory with fasta files to be concatenated", required=True)
    parser.add_argument("-o", "--output", help="Output file name of concatenated sequences", default="concatenated.phylip")
    parser.add_argument("-f", "--format", help="Output sequence format (fasta or phylip)", choices=["phylip","fasta"], default="phylip")
    args = parser.parse_args()

    concat_fams = sorted(glob.glob(f"{args.input}/*.*"))
    all_species = set()
    for fam in concat_fams:
        for r in SeqIO.parse(fam, "fasta"):
            all_species.add(r.id.split("_g")[0])

    all_species = sorted(all_species)
    concatenated_seqs = {sp: "" for sp in all_species}

    for fam in concat_fams:
        recs = list(SeqIO.parse(fam, "fasta"))
        fam_len = len(recs[0].seq)
        seq_by_species = {r.id.split("_g")[0]: str(r.seq) for r in recs}

        for sp in all_species:
            concatenated_seqs[sp] += seq_by_species.get(sp, "-" * fam_len)

    if args.format == "fasta":
        # ------------- write concatenated fasta --------------------------
        with open(args.output, "w") as f:
            for sp in all_species:
                f.write(f">{sp}\n{concatenated_seqs[sp]}\n")
        print(f"[concat] wrote {args.output}")

    if args.format == "phylip":
        # ------------- write concatenated PHYLIP-relaxed -----------------
        records = [SeqRecord(Seq(seq), id=sp, description="") 
                   for sp, seq in concatenated_seqs.items()]
        write_phylip_relaxed(records, args.output)
        print(f"[concat] wrote {args.output} (PHYLIP-relaxed)")

if __name__ == "__main__":
    main()
    exit(0)
