#!/usr/bin/env python3
"""
This script assigns short codes to genomes based on their phyla.

Authors:
    - Lenard Szantho <lenard@drenal.eu>

Version:
    - v0.1 (2025-06-06)
    - v0.2 (2025-06-06) prepare script for duplicate species
"""
import argparse
import csv
import re

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infile", help="Input csv file in the format: genome id,gtdb taxa")
    parser.add_argument("outfile", help="Output file name to save the csv file to. Format will be: genome id,short code")
    args = parser.parse_args()

    shortcodes = {}
    taxa2shortcode = {}
    samespecies = {}
    
    genomes = []
    genome_shortcodes = []
    taxa = []

    with open(args.infile, "r") as inputfh:
        g2t = csv.reader(inputfh)
        for row in g2t:
            genomes.append(row[0])

            if row[1] in taxa2shortcode:
                print("We saw this taxa already {}".format(row[1]))
                if row[1] not in samespecies:
                    # find the first representative
                    i = taxa.index(row[1])
                    genome_shortcodes[i] = "{}s1".format(genome_shortcodes[i])
                samespecies.setdefault(row[1], 1)
                samespecies[row[1]] += 1
                genome_shortcodes.append("{}s{}".format(taxa2shortcode[row[1]], samespecies[row[1]]))
            else:
                shortcode = re.sub(r'fobiota|idibacterota|laeota|bacteria|bacteriota|bacterota|mirabilota|glomotagistota|gistota|glomota|trichota|mirabilota|ficota|trophota|caulota|trophota|phobota|microbiota|mycetota|monadota|dentota|spirota|genetota|chaetota|togota|flexota|oidota|coccota|somatota|vibrionota|cutes|ota', '', row[1].split(";")[1].replace("p__",""))
                shortcode = shortcode.replace("_","").replace("-","")
                if len(shortcode) > 6:
                    shortcode = "{}{}".format(shortcode[:3],shortcode[-3:])

                if re.match("\d", shortcode[-1]):
                    shortcode = shortcode[:-1] + 'n'
                shortcodes.setdefault(shortcode, 0)
                shortcodes[shortcode]+=1

                genome_shortcodes.append("{}{}".format(shortcode, shortcodes[shortcode]))
                taxa2shortcode[row[1]] = genome_shortcodes[-1]

            taxa.append(row[1])

    with open(args.outfile, "w") as outputfh:
        outputfh.write("accession,short,taxa\n")
        for i,v in enumerate(genomes):
            outputfh.write("{},{},{}\n".format(genomes[i], genome_shortcodes[i], taxa[i]))

if __name__ == "__main__":
    main()
