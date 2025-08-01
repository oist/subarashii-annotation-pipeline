#!/usr/bin/env python3
"""
Build families_file.txt for AleRax from the top-300 Phylobayes .treelist or IQtree .ufboot files.

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-23)
"""
import argparse
import sys
import pathlib

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-i", "--input", help="Path to directory containing treelists", required=True, type=pathlib.Path)
    parser.add_argument("-o", "--output", help="Path to output families.txt file", required=True, type=pathlib.Path)
    args = parser.parse_args()

    # try first Phylobayes' treelist
    treelist_paths = list(args.input.glob("*.treelist"))
    if len(treelist_paths) == 0:
        treelist_paths = list(args.input.glob("*.ufboot"))
    if len(treelist_paths) == 0:
        print(f"Neither treelist nor ufboot files found under path {args.input}")
        return 1

    with args.output.open("w") as out:
        out.write("[FAMILIES]\n")
        for p in sorted(treelist_paths):
            fam = p.stem.split(".")[0]          # gf0123
            out.write(f"- {fam}\n")
            out.write(f"gene_tree = {p.resolve()}\n")

    print(f"[prepare_alerax] wrote {args.output}")

if __name__ == "__main__":
    sys.exit(main())

