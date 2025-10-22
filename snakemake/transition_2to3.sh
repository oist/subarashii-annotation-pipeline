#!/bin/bash

usage() {
	printf "\ntransition_2to3.sh script copies a given dataset's\n 2_concatenate_and_filter pipeline's results to the 3_inference\n pipeline's resources directory\n"
	printf "\n Usage: bash transition_2to3.sh <dataset>\n"
}

if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	usage
	exit 0
fi

if [ "$#" -lt 1 ]; then
	echo "[err] This script needs at least one argument."
	usage
	exit 1
fi

if [ "$#" -gt 2 ]; then
	echo "[err] This scripts accepts maximum one argument."
	usage
	exit 2
fi

DATASET=$1
if [ ! -d 2_concatenate_and_filter/results/${DATASET} ]; then
	echo "[err] Directory: 2_concatenate_and_filter/results/${DATASET} doesn't exist, first argument is wrong."
	usage
	exit 3
fi

echo "[log] Creating directory $DATASET for dataset in 3_inference pipeline..."
mkdir -p 3_inference/resources/${DATASET}

cd 3_inference/resources/${DATASET}

for cluster in {mcl,eggnog}; do
	echo "[log] Changing directory to 3_inference/resources/${DATASET}/${cluster}..."
	mkdir ${cluster}
	cd ${cluster}

	echo "[log] Symlinking 2_concatenate_and_filter/results/${DATASET}/${cluster}/concatenated.universality.fasta to 3_inference pipeline's resources directory..."
	ln -s ../../../../2_concatenate_and_filter/results/${DATASET}/${cluster}/concatenated.universality.fa

	echo "[log] Symlinking 2_concatenate_and_filter/results/${DATASET}/${cluster}/concatenated.universality.phylip to to 3_inference pipeline's resources directory..."
	ln -s ../../../../2_concatenate_and_filter/results/${DATASET}/${cluster}/concatenated.universality.phylip

	echo "[log] Creating directory top_families in 3_inference pipeline..."
	mkdir -p top_families

	echo "[log] Changing directory to 3_inference/resources/${DATASET}/${cluster}/top_families..."
	cd top_families

	echo "[log] Symlinking sequences under 2_concatenate_and_filter/results/$DATASET/${cluster}/top_families to 3_inference pipeline's resources directory..."
	for f in ../../../../../2_concatenate_and_filter/results/$DATASET/${cluster}/top_families/*.phylip; do
		ln -s $f
	done
	cd ../..
done

echo "[action required] Please update 3_inference/config/config.yaml with the dataset name: $DATASET to perform the third step of the pipeline on it."

echo "[log] Script ended $(date)."

exit 0
