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
    - v0.1 (2025-10-27)
"""

import argparse
import csv
import statistics
from pathlib import Path
from time import perf_counter
import pyarrow.feather as feather
import pandas as pd



def main():
    t1 = perf_counter()

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--output", help="Output filename", required=True, type=Path)
    parser.add_argument("-m", "--map", help="Mapper csv file between the Eggnog geneIDs and some other IDs to be written to the output file", required=False, type=Path)
    parser.add_argument("--map-col-eggnog", help="Column name corresponding to the values in the Eggnog annotation file")
    parser.add_argument("--map-col-output", help="Column name corresoinding to the output file")
    parser.add_argument("-e", "--eggnog", help="Path to EggNog's annotation file", required=True, type=Path)
    parser.add_argument("-f", "--fullnames", help="GeneID in EggNog annotation is not in the format of speciesID_geneID, the whole name should be used instead of just the speciesID", action="store_true", default=False)
    parser.add_argument("-a", "--all", help="Output all OGs, not just COGs", action="store_true", default=False)
    args = parser.parse_args()

    if args.map:
        if not args.map_col_eggnog or not args.map_col_output:
            print("Missing column definition for the eggnog column")
            exit(2)

    # stores COG counts per genomes
    # genome1 : { COG0004: 1, COG0006: 4, ... }
    genome2cogcounts = dict()
    cogs = set()

    print(f"Processing EggNog annotation file {args.eggnog} ...")
    with open(args.eggnog) as csvfile:
        annotations = csv.reader(csvfile, delimiter='\t')
        for rows in annotations:
            if rows[0][0] == "#":
               continue
            gene = rows[0].strip()
            genome = gene
            if not args.fullnames:
                genome = gene.split("_")[0]
            cog = rows[4].split("@")[0].strip()

            if not args.all:
                # filter out non COGs
                if cog[:3] == "COG":
                    genome2cogcounts.setdefault(genome, {})
                    genome2cogcounts[genome].setdefault(cog, 0)
                    genome2cogcounts[genome][cog] += 1
                    cogs.add(cog)

    # Create dataframe
    genome2cogcounts_df = pd.DataFrame.from_dict(genome2cogcounts, orient="index")
    # Sort COGs by their name
    genome2cogcounts_df = genome2cogcounts_df[sorted(genome2cogcounts_df.columns.values)]
   
    # change names if mapping is provided
    if args.map:
        print(f"Loading mapping file {args.map} ...")
        old_columns_list = sorted(genome2cogcounts_df.columns.values)
        mapping = pd.read_csv(args.map, skipinitialspace=True)
        genome2cogcounts_df = pd.merge(left=genome2cogcounts_df, right=mapping, left_index=True, right_on=args.map_col_eggnog, how="left")
        genome2cogcounts_df = genome2cogcounts_df[ sorted(mapping.columns.values) + old_columns_list ]
        # not nice, but works:
        #genome2cogcounts_df[args.map_col_output] = genome2cogcounts_df.index
        #genome2cogcounts_df[args.map_col_output] = genome2cogcounts_df[args.map_col_output].replace(to_replace={ k:v for k, v in mapping.loc[:,[args.map_col_eggnog,args.map_col_output]].values})

    print(f"Writing presence-absence table to {args.output} ...")
    index_label = "genome"
    if args.map:
        index_label=args.map_col_eggnog

    genome2cogcounts_df.to_csv(args.output, header=True, index=False)
    feather.write_feather(genome2cogcounts_df, f"{args.output}.feather")


if __name__ == "__main__":
    main()
    exit(0)

