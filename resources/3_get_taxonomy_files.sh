#!/bin/bash

wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/bac120_taxonomy_r202.tsv.gz

gzip -d bac120_taxonomy_r202.tsv.gz

REGEXP=`tr '\n' '|' < list_of_genomes_accession.txt  | sed 's/|$//'`

grep -E "$REGEXP" bac120_taxonomy_r202.tsv | sed -E 's/\s+/,/g' > genome2taxa.csv
# adding manually one of the genomes that is not present in GTDB r202
echo "GCF_004343465.1,d__Bacteria;p__Firmicutes_A;c__Thermoanaerobacteria;o__Thermoanaerobacterales;f__Thermoanaerobacteraceae;g__Thermoanaerobacterium;s__Thermoanaerobacterium thermosaccharolyticum" >> genome2taxa.csv



