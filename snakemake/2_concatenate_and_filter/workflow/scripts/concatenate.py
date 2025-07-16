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
import sys, subprocess, statistics, pathlib, shutil
from collections import defaultdict
import pandas as pd
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
    parser.add_argument("-i", "--input", help="Directory with fasta files or list of fasta files to be concatenated", required=True, nargs="+")
    parser.add_argument("-o", "--output", help="Output file name of concatenated sequences", default="concatenated")
    parser.add_argument("-f", "--format", help="Sequence format (fasta or phylip)", choices=["phylip","fasta"], default="phylip")
    parser.add_argument("-t", "--threshold", help="Universality threshold", default=0.8, type=float)
    parser.add_argument("-m", "--metadata", help="Metadata tsv or csv file containing a per-sequence data")
    parser.add_argument("-c", "--columnname", help="Name of column in metadata file for threshold", default="universality")
    args = parser.parse_args()

    # --- load scores --------------------------------------------------
    df = pd.read_csv(SCORE_TSV, sep="\t")
    df["universality"] = df["n_species"] / df["n_species"].max()

    # ---------- 1. Concatenation set ---------------------------------
    concat_fams = df[df["universality"] >= 0.80]["family"].tolist()
    print(f"[concat] {len(concat_fams)} families ≥80 % universality")

    # first pass: discover universe of species
    all_species = set()
    for fam in concat_fams:
        trim = ALIGN_DIR / f"{fam}.trim"
        if not trim.exists():
            continue
        for r in SeqIO.parse(trim, "fasta"):
            all_species.add(r.id.split("_g")[0])

    all_species = sorted(all_species)
    concatenated_seqs = {sp: "" for sp in all_species}

    for fam in concat_fams:
        trim = ALIGN_DIR / f"{fam}.trim"
        if not trim.exists():
            continue

        # copy .trim once
        dst = CONCAT_DIR / trim.name
        if not dst.exists():
            shutil.copy(trim, dst)

        recs = list(SeqIO.parse(trim, "fasta"))
        fam_len = len(recs[0].seq)
        seq_by_species = {r.id.split("_g")[0]: str(r.seq) for r in recs}

        for sp in all_species:
            concatenated_seqs[sp] += seq_by_species.get(sp, "-" * fam_len)

    # ------------- write concatenated fasta --------------------------
    concat_faa = CONCAT_DIR / "concatenated.faa"
    with concat_faa.open("w") as f:
        for sp in all_species:
            f.write(f">{sp}\n{concatenated_seqs[sp]}\n")
    print(f"[concat] wrote {concat_faa}")

    # ------------- write concatenated PHYLIP-relaxed -----------------
    concat_phy = CONCAT_DIR / "concatenated.phy"
    records = [SeqRecord(Seq(seq), id=sp, description="") 
               for sp, seq in concatenated_seqs.items()]
    write_phylip_relaxed(records, concat_phy)
    print(f"[concat] wrote {concat_phy} (PHYLIP-relaxed)")

if __name__ == "__main__":
    main()
    exit(0)
