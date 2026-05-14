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
    - v0.7 (2026-03-19)
        - feather export out of memory for 150k x 207k matrix, trying to make more memory-efficient
    - v0.6 (2026-03-16)
        - feather export restored and made into an option
    - v0.5 (2026-03-15)
        - bugfix: one extra comma at the end of lines removed
    - v0.4 (2026-02-26)
        - generalizing for KEGG KO and EC category columns as well
    - v0.3 (2026-02-13)
        - bugfix: count all |root COGs, not just the first one
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
from memory_profiler import profile

def main():
    t1 = perf_counter()

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-e", "--eggnog", help="Path to EggNog's annotation file", required=True, type=Path)
    parser.add_argument("-o", "--output", help="Output filename base", required=True, type=Path)
    parser.add_argument("-m", "--map", help="Metadata mapping csv file with 2 or more columns to be included in the genomes2cogcounts file.", required=False, type=Path)
    parser.add_argument("--map-col-eggnog", help="Column name in the mapping csv file corresponding to the genome column in the eggnog annotation file.")
    parser.add_argument("-a", "--all", help="Output all OGs, not just COGs. Has effect only if '-c' is set to 'cog'", action="store_true", default=False)
    parser.add_argument("-c", "--column", help="Column to count. By default the 'cog', i.e. counting the eggnog_OG column.", default="cog", choices=['cog', 'ko', 'ec'])
    parser.add_argument("-f", "--feather", help="Output as feather format as well", action="store_true", default=False)
    parser.add_argument("-r", "--reduce", help="Reduce the genomes to a list of accession IDs or genome IDs (supply a file with a list separated by newline chars)")
    args = parser.parse_args()

    columns_mapping = {
            "cog": 4,
            "ko": 11,
            "ec": 10
    }

    if args.map:
        if not args.map_col_eggnog:
            print("Error: Missing mapping csv file column definition corresponding to the genome column in the eggnog annotation file")
            exit(2)

    # stores COG counts per genomes
    # { "genome1" : { "COG0004": 1, "COG0006": 4, ... } }
    genome2cogcounts = dict()
    cogs = set()
    # stores COG association of a gene (same as input basically)
    # { "gene1": ["COG0004", "COGXXXX", ...] }
    gene2cogs = dict()
    # stores genes associated to a given COG
    # { "COG0043": ["gene1", "gene2", "gene3", ...] }
    cog2genes = dict()
    # genomes associated with a cog
    cog2nrgenomes = dict()


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
            cog_list = []
            # e.g.: COG0457@1|root,COG0823@1|root,COG2885@1|root,COG0457@2|Bacteria,COG0823@2|Bacteria,COG2885@2|Bacteria,4NE6G@976|Bacteroidetes,2FPQX@200643|Bacteroidia
            for c in rows[columns_mapping[args.column]].split(","):
                if args.column == "cog" and "|root" in c:
                    # e.g. COG0457@1|root
                    if args.all or c[:3] == "COG":
                        cog_list.append(c.split("@")[0].strip())
                if args.column == "ko":
                    #  ko:K02052,ko:K11072,ko:K11076
                    cog_list.append(c[3:].strip())
                if args.column == "ec":
                    cog_list.append(c.strip())

            genome2cogcounts.setdefault(genome, {})
            for c in cog_list:
                genome2cogcounts[genome].setdefault(c, 0)
                genome2cogcounts[genome][c] += 1

                cog2nrgenomes.setdefault(c,0)
                cog2nrgenomes[c] += 1

                cogs.add(c)

                cog2genes.setdefault(c, [])
                cog2genes[c].append(gene)
            gene2cogs[gene] = cog_list

    # change names if mapping is provided
    if args.map:
        print(f"Loading mapping file {args.map} ...")
        # { "genomeID" : {"accession" : "RS_GCF_01XXXXX", "taxa": "d__Bacteria;p__;c__;p__;f__;g__;s__" }, "genomeID2": {}, ... }
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

    if args.map and args.reduce:
        accessions_to_keep = []
        with open(args.reduce, "r") as inputfh:
            for l in inputfh:
                if len(l.strip()) > 0:
                    accessions_to_keep.append(l.strip())
        print(f"{len(accessions_to_keep)} number of accessions should be retained.")
        genome_ids_to_delete = []
        for genome in genome2cogcounts:
            if mapping[genome]["accession"] not in accessions_to_keep:
                genome_ids_to_delete.append(genome)
        print(f"{len(genome_ids_to_delete)} genomes will be removed from the count table")
        for g in genome_ids_to_delete:
            if g in genome2cogcounts:
                del genome2cogcounts[g]

    # check how often we have singleton COGs
    cog_to_delete = []
    for cog in cog2nrgenomes:
        if cog2nrgenomes[cog] == 1:
            #print(f"{cog} has only one genome in it")
            cog_to_delete.append(cog)

    print(f"{len(cog_to_delete)} categories are singleton, to be deleted")

    for cog in cog_to_delete:
        for genome in genome2cogcounts:
            if cog in genome2cogcounts[genome]:
                del genome2cogcounts[genome][cog]
        cogs.discard(cog)
        

    print(f"Writing count table to {args.output} ...")
    index_label = "genome"

    with open(f"{args.output}_genome2{args.column}counts.csv", "w") as outputfh:
        if args.map:
            outputfh.write(f"{index_label},accession,taxa,")
        else:
            outputfh.write(f"{index_label},")
        outputfh.write(",".join(sorted(cogs)))
        outputfh.write("\n")
        for genome in genome2cogcounts:
            outputfh.write(f"{genome}")
            if args.map:
                outputfh.write(f",{mapping[genome]['accession']},{mapping[genome]['taxa']}")
            for cog in sorted(cogs):
                if cog in genome2cogcounts[genome]:
                   outputfh.write(f",{genome2cogcounts[genome][cog]}")
                else:
                    outputfh.write(f",")
            outputfh.write("\n")

    with open(f"{args.output}_{args.column}2genes.list", "w") as outputfh:
        for cog in sorted(cog2genes.keys()):
            outputfh.write("{}\t{}\n".format(cog, ",".join(cog2genes[cog])))

    if args.feather:
        import pyarrow.feather as feather
        import pandas as pd
        import numpy as np

        # trying to use intc to reduce space, but apparently pandas just ignores it and uses doubles
        df = pd.DataFrame(columns=sorted(cogs), index=list(genome2cogcounts.keys()), dtype=np.intc)
        for genome in genome2cogcounts:
            df.loc[genome] = pd.Series(genome2cogcounts[genome])
        # this duplicates the DataFrame in the memory, that would lead to OOF-kill for large datasets
        #df = df.fillna(0)

        feather.write_feather(df, f"{args.output}_genome2{args.column}counts.feather")


if __name__ == "__main__":
    main()
    exit(0)

