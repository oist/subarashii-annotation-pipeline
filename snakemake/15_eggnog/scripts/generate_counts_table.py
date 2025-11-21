#!/usr/bin/env python3
"""
Generate COG presence-absence table based on EggNog's annotation file

TODO:
    - new output also: COG\tgene1,gene2,geneN
    - count table: list members instead of just the number
    - sqlite for the tasks of: all sequence in a genome, all fasta in a COG, etc:
        - genome, name, gtdb metadata, protein coding sequences, COG
        - would it be faster than pandas, etc?
        - how can one of these be hosted on Deigo?
        - how many genes are COGs?
    - embedding: model with 1000 (less) dimensions, but much faster
Authors:
    - Lenard Szantho <lenard@drenal.eu>

Version:
    - v0.2 (2025-11-18) 
        - dropping pandas and pyarrow, due to RAM usage, using just dictionaries, 
        - exporting cog2genes file with the structure: "COG_ID    gene1,gene2,gene3,...,geneN"
    - v0.1 (2025-10-27)
        - feather file export with pyarrow
        - genome2cogcounts file
"""

import argparse
import csv
import statistics
from pathlib import Path
from time import perf_counter
#import pyarrow.feather as feather
#import pandas as pd



def main():
    t1 = perf_counter()

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-e", "--eggnog", help="Path to EggNog's annotation file", required=True, type=Path)
    parser.add_argument("-o", "--output", help="Output filename base", required=True, type=Path)
    parser.add_argument("-m", "--map", help="Metadata mapping csv file with 2 or more columns to be included in the genomes2cogcounts file.", required=False, type=Path)
    parser.add_argument("--map-col-eggnog", help="Column name in the mapping csv file corresponding to the genome column in the eggnog annotation file.")
    parser.add_argument("-a", "--all", help="Output all OGs, not just COGs", action="store_true", default=False)
    args = parser.parse_args()

    if args.map:
        if not args.map_col_eggnog:
            print("Error: Missing mapping csv file column definition corresponding to the genome column in the eggnog annotation file")
            exit(2)

    # stores COG counts per genomes
    # { "genome1" : { "COG0004": 1, "COG0006": 4, ... } }
    genome2cogcounts = dict()
    cogs = set()
    # stores COG association of a gene (same as input basically)
    # { "gene1": "COG0004" }
    gene2cogs = dict()
    # sotores genes associated to a given COG
    # { "COG0043": ["gene1", "gene2", "gene3", ...] }
    cog2genes = dict()

    print(f"Processing EggNog annotation file {args.eggnog} ...")
    with open(args.eggnog) as csvfile:
        annotations = csv.reader(csvfile, delimiter='\t')
        for rows in annotations:
            if rows[0][0] == "#":
                continue
            # gene: GenomeShortCode_geneID
            # e.g. Desulfo10_g1909
            gene = rows[0].strip()
            genome = gene.split("_")[0]
            cog = rows[4].split("@")[0].strip()
            
            if not args.all and not cog[:3] == "COG":
                continue

            genome2cogcounts.setdefault(genome, {})
            genome2cogcounts[genome].setdefault(cog, 0)
            genome2cogcounts[genome][cog] += 1
            cogs.add(cog)

            cog2genes.setdefault(cog, [])
            cog2genes[cog].append(gene)
            gene2cogs[gene] = cog

    # change names if mapping is provided
    if args.map:
        print(f"Loading mapping file {args.map} ...")
        # { "genomeID" : {"accession" : "RS_GCF_01XXXXX", "taxa": "d__Bacteria;p__;c__;p__;f__;g__;s__"
        mapping = {}
        mapping_columns = []
        index_col_id = 0
        accession_col_id = 0
        taxa_col_id = 0
        with open(args.map, "r") as csvfile:
            mappings = csv.reader(csvfile)
            firstrow = True
            for row in mappings:
                if firstrow:
                    firstrow = False
                    mapping_columns = row
                    if not args.map_col_eggnog in mapping_columns:
                        print(f"Error: {args.map_col_eggnog} not found in the mapping {args.map} file.")
                        exit(5)
                    index_col_id = mapping_columns.index(args.map_col_eggnog)
                    accession_col_id = mapping_columns.index("accession")
                    taxa_col_id = mapping_columns.index("taxa")

                mapping.setdefault(row[index_col_id], {})
                mapping[row[index_col_id]]["accession"] = row[accession_col_id]
                mapping[row[index_col_id]]["taxa"] = row[taxa_col_id]

    print(f"Writing count table to {args.output} ...")
    index_label = "genome"

    with open(f"{args.output}_genome2cogcounts.csv", "w") as outputfh:
        outputfh.write(f"{index_label},accession,taxa,")
        outputfh.write(",".join(sorted(cogs)))
        outputfh.write("\n")
        for genome in genome2cogcounts:
            outputfh.write(f"{genome},")
            if args.map:
                outputfh.write(f"{mapping[genome]['accession']},{mapping[genome]['taxa']},")
            for cog in sorted(cogs):
                if cog in genome2cogcounts[genome]:
                    outputfh.write(f"{genome2cogcounts[genome][cog]},")
                else:
                    outputfh.write(",")
            outputfh.write("\n")

    with open(f"{args.output}_cog2genes.list", "w") as outputfh:
        for cog in sorted(cog2genes.keys()):
            outputfh.write("{}\t{}\n".format(cog, ",".join(cog2genes[cog])))


if __name__ == "__main__":
    main()
    exit(0)

