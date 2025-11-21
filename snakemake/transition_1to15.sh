#!/bin/bash

usage() {
	printf "\ntransition_1to15.sh script copies a given dataset's\n 1_cluster pipeline's selected genomes to the 15_eggnog\n pipeline's resources directory\n"
	printf "In order for a more efficient eggnog annotation process\n"
	printf "\n Usage: bash transition_1to15.sh <dataset> <parts>\n"
}

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	usage
	exit 0
fi

if [ "$#" -lt 2 ]; then
	echo "[err] This script needs at least 2 arguments."
	usage
	exit 1
fi

if [ "$#" -gt 2 ]; then
	echo "[err] This scripts accepts maximum 2 arguments."
	usage
	exit 2
fi

DATASET=$1
if [ ! -d 1_cluster/results/${DATASET} ]; then
	echo "[err] Directory: 1_cluster/results/${DATASET} doesn't exist, first argument is wrong."
	usage
	exit 3
fi

if [ ! -f 1_cluster/results/${DATASET}/all_genomes_protein.faa ]; then
	echo "[err] Directory: 1_cluster/results/${DATASET}/all_genomes_protein.faa doesn't exist."
	usage
	exit 4
fi

PARTS=$2
if [ "$PARTS" -lt "1" ]; then
	echo "[err] Parts has to be 1 or more"
	usage
	exit 5
fi

echo "[log] Creating directory $DATASET for dataset in 15_eggnog pipeline ..."
mkdir -p 15_eggnog/resources/${DATASET}

echo "[log] Changing directory to 15_eggnog/resources/${DATASET}"
cd 15_eggnog/resources/${DATASET}

echo "[log] Counting number of sequences in all_genomes_protein.faa . May take several seconds, please be patient."
#NUMBER_OF_SEQS=`grep ">" ../../../1_cluster/results/${DATASET}/all_genomes_protein.faa | wc -l`
NUMBER_OF_LINES=`wc -l ../../../1_cluster/results/${DATASET}/all_genomes_protein.faa | sed -E 's/\s+/ /' | cut -d " " -f 1`
# we assume that each sequence is in one line, thus half of the lines in the file is the number of sequences
BATCH_SIZE=`awk "function ceil(x){return int(x)+(x>int(x))} BEGIN { print 2 * ceil( ( $NUMBER_OF_LINES / 2 ) / $PARTS) }"`

echo "[log] Breaking all_genomes_protein.faa into $PARTS parts of size $((BATCH_SIZE / 2)) sequences ..."
for i in `seq 0 $PARTS`; do
	sed -n "$((i*BATCH_SIZE+1)),+$((BATCH_SIZE-1))p" ../../../1_cluster/results/${DATASET}/all_genomes_protein.faa > all_genomes_protein.part${i}.faa
done

echo "[log] Symlinking genome2abbrev.csv ..."
ln -s ../../../1_cluster/results/${DATASET}/genome2abbrev.csv

echo "[log] Symlinking eggnog_diamond_db ..."
cd ..
mkdir eggnog_diamond_db
for f in ../../../1_cluster/resources/eggnog_diamond_db/*; do
	ln -s $f
done

echo "[action required] Please update 15_eggnog/config/config.yaml with the new dataset name: $NEW_DATASET to perform the splitted eggnog annotation step of the pipeline on it."

echo "[log] Script ended $(date)."

exit 0
