#!/usr/bin/env python3
"""
Create COG clusters based on EggNog annotation file

Also writes four reports in the same directory:
   cog2genes.tsv
   genes2cog.tsv
   family_type_counts.tsv -- Singleton, Duplets, Triplets and Huge counts
   family_details.tsv -- number of sequences and sequence length stats and if it has paralogs or not

Authors:
    - Lenard Szantho <lenard@drenal.eu>

Version:
    - v0.2 (2026-05-14): handling multiple root-level COGs
    - v0.1 (2025-10-07)
"""

import argparse
import csv
import statistics
from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO
from collections import Counter
from time import perf_counter


def has_paralogs(genes):
    species_seen = Counter(gene.split("_g")[0] for gene in genes)
    return any(cnt > 1 for cnt in species_seen.values())

def main():
    t1 = perf_counter()

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", help="Output base directory.", required=True)
    parser.add_argument("-f", "--fasta", help="Directory with the per genome fasta files", required=True)
    parser.add_argument("-e", "--eggnog", help="Path to EggNog's annotation file", required=True)
    parser.add_argument("-c", "--cutoff", help="Cutoff number for huge families. By default: 5 times the number of genomes.", type=int)
    args = parser.parse_args()

    out_root = Path(args.output)
    out_root.mkdir(parents=True, exist_ok=True)

    genome2genes2cog = dict()
    cog2genes = dict()
    cog2seqlengths = dict()
    print(f"Huge cutoff set to {args.cutoff}")
    print(f"Processing EggNog annotation file {args.eggnog} ...")
    with open(args.eggnog) as csvfile:
        annotations = csv.reader(csvfile, delimiter='\t')

        # Detect column indices from the header line (last "#" line before data).
        # v2: eggNOG_OGs holds "COGxxxx@1|root,..." — filter for |root entries.
        # v3: eggNOG_OGs holds domain names; COG_category holds the COG ID directly.
        og_col = 4        # eggNOG_OGs, v2 default
        cog_cat_col = None  # COG_category (v3: full COG ID; v2: single-letter code)
        for rows in annotations:
            if rows[0].startswith("#"):
                header = [h.strip() for h in ([rows[0].lstrip("#").strip()] + rows[1:])]
                for i, h in enumerate(header):
                    if h == "eggNOG_OGs":
                        og_col = i
                    elif h == "COG_category":
                        cog_cat_col = i
                continue
            break

        print(f"Column indices — eggNOG_OGs: {og_col}, COG_category: {cog_cat_col}")

        csvfile.seek(0)
        annotations = csv.reader(csvfile, delimiter='\t')
        for rows in annotations:
            if rows[0][0] == "#":
                continue
            gene = rows[0].strip()
            genome = gene.split("_g")[0]

            cog_list = []
            # v2: COG IDs appear as "COGxxxx@1|root" inside the eggNOG_OGs column.
            # e.g.: COG0457@1|root,COG0823@1|root,COG2885@1|root,COG0457@2|Bacteria
            for c in rows[og_col].split(","):
                if "|root" in c:
                    token = c.split("@")[0].strip()
                    if token.startswith("COG"):
                        cog_list.append(token)

            # v3 fallback: COG_category column now contains the full COG ID directly.
            if not cog_list and cog_cat_col is not None and cog_cat_col < len(rows):
                cog_cat = rows[cog_cat_col].strip()
                if cog_cat and cog_cat != "-" and cog_cat.startswith("COG"):
                    cog_list.append(cog_cat)

            genome2genes2cog.setdefault(genome, {})
            genome2genes2cog[genome].setdefault(gene, [])
            for c in cog_list:
                genome2genes2cog[genome][gene].append(c)
                cog2genes.setdefault(c, [])
                cog2genes[c].append(gene)

    print(f"Sorting {len(genome2genes2cog)} genomes into COG files under {out_root} ...")
    for genome in tqdm(genome2genes2cog):
        genome_file = {rec.id: rec for rec in SeqIO.parse(f"{args.fasta}/{genome}.faa", "fasta")}
        for gene in genome2genes2cog[genome]:
            for c in genome2genes2cog[genome][gene]:
                cog2seqlengths.setdefault(c, []).append(len(genome_file[gene].seq))
                with open(out_root / f"{c}.faa", "a") as outputfh:
                    outputfh.write(f">{gene}\n{genome_file[gene].seq}\n")

    print("Order families into subdirectories based on the number of genomes")
    type_counter = dict()
    for t in ("Singleton", "Duplets", "Triplets", "Normal", "Huge"):
        ( out_root / t ).mkdir(exist_ok=True)
        type_counter[t] = 0

    cutoff_size = 300
    if not args.cutoff:
        cutoff_size = 5 * len(genome2genes2cog)
    else:
        cutoff_size = args.cutoff

    for cog in cog2genes:
        size = len(cog2genes[cog])
        if size == 1:
            type_counter["Singleton"] += 1
            ( out_root / 'Singleton' / f"{cog}.faa" ).symlink_to(f"../{cog}.faa")
        elif size == 2:
            type_counter["Duplets"] += 1
            ( out_root / 'Duplets' / f"{cog}.faa" ).symlink_to(f"../{cog}.faa")
        elif size == 3:
            type_counter["Triplets"] += 1
            ( out_root / 'Triplets' / f"{cog}.faa" ).symlink_to(f"../{cog}.faa")
        elif size >= cutoff_size:
            type_counter["Huge"] += 1
            ( out_root / 'Huge' / f"{cog}.faa" ).symlink_to(f"../{cog}.faa")
        else:
            type_counter["Normal"] += 1
            ( out_root / 'Normal' / f"{cog}.faa" ).symlink_to(f"../{cog}.faa")

    # write logs
    with open(out_root / 'genes2cog.tsv', 'w') as outputfh:
        for genome in genome2genes2cog:
            for gene in genome2genes2cog[genome]:
                cogs = ";".join(genome2genes2cog[genome][gene])
                outputfh.write(f"{gene}\t{cogs}\n")

    with open(out_root / 'cog2genes.tsv', 'w') as outputfh:
        for cog in cog2genes:
            outputfh.write(f"{cog}\t")
            outputfh.write(",".join(cog2genes[cog]))
            outputfh.write("\n")
    
    with open(out_root / 'family_type_counts.tsv', 'w') as outputfh:
        #outputfh.write("Name\tThreshold\tSize\n")
        for t in ("Singleton", "Duplets", "Triplets", "Normal", "Huge"):
            outputfh.write(f"{t}\t{type_counter[t]}\n")

    with open(out_root / 'family_details.tsv', 'w') as outputfh:
        outputfh.write("family\tsize\tlongest\tshortest\tmedian\tparalog_in_species\n")
        for cog in cog2genes:
            size = len(cog2genes[cog])
            max_seqlen = max(cog2seqlengths[cog])
            min_seqlen = min(cog2seqlengths[cog])
            median_seqlen = statistics.median(cog2seqlengths[cog])
            paralogs = "Yes" if has_paralogs(cog2genes[cog]) else "No"
            outputfh.write(f"{cog}\t{size}\t{max_seqlen}\t{min_seqlen}\t{median_seqlen}\t{paralogs}\n")

    normal_count = type_counter.get("Normal", 0)
    if normal_count == 0:
        print(
            f"[create_cog_clusters] ERROR: 0 Normal-sized COG families produced.\n"
            f"  eggNOG_OGs column used: index {og_col}\n"
            f"  This usually means the column index is wrong for your EggNOG-mapper version.\n"
            f"  Check the header of {args.eggnog} and ensure eggNOG_OGs is at index {og_col}.",
            flush=True,
        )
        raise SystemExit(1)

    t2 = perf_counter()
    print(f"[create_cog_clusters] done. Elapsed time: {t2-t1:.2f} secs. Outputs in: {out_root}")


if __name__ == "__main__":
    main()
    exit(0)

