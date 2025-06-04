#!/bin/bash

echo "This script is not needed, as we'll use the aa sequences."
exit 0

mkdir genomes1007

cd genomes1007

for f in `cat ../list_of_genomes_accession.txt | awk -F "" '{print $4$5$6"/"$8$9$10"/"$11$12$13"/"$14$15$16"/"$4$5$6"_"$8$9$10$11$12$13$14$15$16"."$18"_genomic.fna.gz"}'`; do
	if [ -f ../gtdb_genomes_reps_r202/$f ]; then
		cp ../gtdb_genomes_reps_r202/$f ./
	else
		echo "$f not found"
	fi
done
