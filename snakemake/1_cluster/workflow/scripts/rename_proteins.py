#!/usr/bin/env python3
"""
Rename proteins in FASTA files for the subarashii pipeline.

Outputs written to --out
    <ShortCode>.faa          renamed FASTA (headers: >ShortCode_gN)
    genome2abbrev.tsv        accession,short[,taxa]
    geneid_mapping.txt       new_gene_id old_full_header
    mapping_autogen.csv      only when --auto is used

Auhtors:
    - Adrian A. Davin
"""

import csv
import pathlib
import argparse
import sys
from Bio import SeqIO
from tqdm import tqdm
import glob

# ---------- argument parsing ----------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("-m", "--map", type=pathlib.Path,
                   help="2 or 3 column CSV/TSV: original_faa,ShortCode[,taxa]")
    p.add_argument("-a", "--auto", action="store_true",
                   help="ignore --map and generate s0, s1, ...")
    p.add_argument("-i", "--in", dest="infolder", required=True, type=pathlib.Path)
    p.add_argument("-o", "--out", dest="outfolder", required=True, type=pathlib.Path)
    p.add_argument("-t", "--translationoutfile", help="Filename in which the translation between renamed and original protein IDs will be written", required=True)
    args = p.parse_args()
    if not args.auto and args.map is None:
        p.error("provide --map or add --auto")
    return args


# ---------- mapping helpers ----------
def generate_auto_mapping(folder):
    faa_files = sorted(folder.glob("*.faa"))
    return [(f.name, f"s{i}", "") for i, f in enumerate(faa_files)]


def read_mapping(path):
    rows = []
    with path.open() as fh:
        dialect = csv.Sniffer().sniff(fh.readline())
        fh.seek(0)
        for cols in csv.reader(fh, dialect):
            if len(cols) == 2:
                rows.append((cols[0], cols[1], ""))
            elif len(cols) >= 3:
                rows.append((cols[0], cols[1], cols[2]))
            else:
                raise ValueError("Bad row in mapping file: " + str(cols))
    return rows


def write_csv(rows, dest_path):
    with dest_path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter=",")
        for acc, short, taxa in rows:
            if taxa:
                writer.writerow([acc, short, taxa])
            else:
                writer.writerow([acc, short])


# ---------- renaming core ----------
def rename_and_map(in_faa, out_faa, prefix, gene_map_fh):
    total = sum(1 for _ in SeqIO.parse(in_faa, "fasta"))
    records = []
    for i, rec in enumerate(
            tqdm(SeqIO.parse(in_faa, "fasta"),
                 total=total,
                 desc=in_faa,
                 unit="seq"), 1):

        old_header = rec.description       # full original header (no leading >)
        new_id = f"{prefix}_g{i}"

        rec.id = new_id
        rec.description = ""               # header will be >new_id only

        # TODO csv or tsv?
        # TODO include file name (accession ID) as well?
        gene_map_fh.write(f"{new_id}\t{old_header}\n")
        records.append(rec)

    SeqIO.write(records, out_faa, "fasta")


# ---------- main ----------
def main():
    args = parse_args()
    args.outfolder.mkdir(parents=True, exist_ok=True)

    # build mapping list
    if args.auto:
        mapping = generate_auto_mapping(args.infolder)
        write_tsv(mapping, args.outfolder / "mapping_autogen.csv")
    else:
        mapping = read_mapping(args.map)

    missing_files = []
    with open(args.translationoutfile, "w", newline="") as gm:

        g2a_rows = []
        for original, short, taxa in tqdm(mapping,
                                          desc="FASTA files",
                                          unit="file"):

            in_faa = glob.glob(r"{}/{}*.fa*".format(args.infolder, original))
            if len(in_faa) < 1:
                missing_files.append(args.infolder / original)
                continue
            if len(in_faa) > 1:
                print("Accession column {} in mapping file doesn't uniquely determine sequences filename. Exiting...".format(args.infolder / original))
                exit(2)

            out_faa = args.outfolder / f"{short}.faa"
            rename_and_map(in_faa[0], out_faa, short, gm)

            accession = original.rsplit(".", 1)[0]   # drop .faa
            g2a_rows.append((accession, short, taxa))


    if len(missing_files) > 1:
        # if only one, probably it's the header line, if more, it's a problem
        print("Error: Check your data. The following genome files that were in the mapping file could not be found:")
        for f in missing_files:
            print(f)
        return 1


if __name__ == "__main__":
    sys.exit(main())

