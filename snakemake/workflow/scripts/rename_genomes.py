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
    parser.add_argument("-d", "--domain", help="Prefix species IDs with domain's first letter", action="store_true", default=False)
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
                shortcode = re.sub(r'archaeota|archaeum|fobiota|idibacterota|laeota|bacteria|bacteriota|bacterota|mirabilota|glomotagistota|gistota|glomota|trichota|mirabilota|ficota|trophota|caulota|trophota|phobota|microbiota|mycetota|monadota|dentota|genetota|chaetota|togota|flexota|oidota|coccota|somatota|vibrionota|cutes|microbia|ota', '', row[1].split(";")[1].replace("p__",""))

                domain = row[1].split(";")[0].replace("d__","")[0]
                shortcode = shortcode.replace("_","").replace("-","")
                shortcode = shortcode.replace("Thermo","The").replace("stone","st")
                #if len(shortcode) > 6:
                #    #shortcode = "{}{}".format(shortcode[:3],shortcode[-3:])
                #    shortcode = "{}".format(shortcode[:6])

                # we don't need to prefix with non-number,
                # if it gets prefixed by domain's first letter
                if not args.domain and re.match("\d", shortcode[0]):
                    shortcode = 'S' + shortcode

                if re.match("\d", shortcode[-1]):
                    if re.match("\d", shortcode[-3]):
                        shortcode = shortcode[:-2] + 'n'
                    else:
                        shortcode = shortcode[:-2]
                shortcodes.setdefault(shortcode, 0)
                shortcodes[shortcode]+=1

                if args.domain:
                    genome_shortcodes.append("{}{}{}".format(domain, shortcode, shortcodes[shortcode]))
                else:
                    genome_shortcodes.append("{}{}".format(shortcode, shortcodes[shortcode]))
                taxa2shortcode[row[1]] = genome_shortcodes[-1]

            taxa.append(row[1])

    with open(args.outfile, "w") as outputfh:
        outputfh.write("accession,short,taxa\n")
        for i,v in enumerate(genomes):
            outputfh.write("{},{},{}\n".format(genomes[i], genome_shortcodes[i], taxa[i]))

if __name__ == "__main__":
    main()
