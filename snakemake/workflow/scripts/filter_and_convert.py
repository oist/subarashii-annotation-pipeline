#!/usr/bin/env python3
"""
Read DIAMOND hits (tab-separated), keep only matches that pass the
e-value and coverage cut-offs and write the result in MCL .abc 
format (q s weight).

Authors:
    - Adrian A. Davin

Version:
    - v0.1 (2025-07-03)
"""
import csv
import argparse

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-e", "--evalue", help="e-value cut-off", required=True)
    parser.add_argument("-c", "--coverage", help="Coverage cut-off (in percentage)", required=True)
    parser.add_argument("-i", "--input", help="diamond's all against all output tsv file", required=True)
    parser.add_argument("-o", "--output", help="Output file name for MCL-compatible abc file", required=True)
    args = parser.parse_args()

    kept = 0
    total = 0

    with args.input.open() as fi, args.output.open("w") as fo:
        reader = csv.reader(fi, delimiter="\t")
        for q, s, pid, bits, e, qcov, scov in reader:
            total += 1
            if float(e) > args.evalue:
                continue
            if min(float(qcov), float(scov)) < args.coverage:
                continue
            fo.write(f"{q} {s} {bits}\n")
            kept += 1

    print(f"[filter_and_convert] kept {kept}/{total} hits "
          f"(E <= {args.evalue}, cov >= {args.coverage}%)")

if __name__ == "__main__":
    main()
    exit(0)

