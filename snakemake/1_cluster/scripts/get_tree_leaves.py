#!/usr/bin/env python3
"""
This scripts outputs the list of leaf names from a tree file.

Authors:
    Lenard Szantho <lenard@drenal.eu>

Version:
    - v0.1 (2025-08-29)
"""
import sys
import argparse
from newick import read

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("treefile", help="Newick tree file")
    parser.add_argument("-o", "--output", help="Output file (if not set, stdout)")
    args = parser.parse_args()

    t = read(args.treefile)[0]
    leaves = t.get_leaf_names()

    if args.output:
        with open(args.output, "w") as outputfh:
            for l in leaves:
                outputfh.write(f"{l}\n")
    else:
        for l in leaves:
            print(l)


if __name__ == "__main__":
    sys.exit(main())
