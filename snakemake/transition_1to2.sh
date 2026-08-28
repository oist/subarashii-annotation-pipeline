#!/bin/bash

usage() {
	printf "\ntransition_1to2.sh script copies a given dataset's\n 1_cluster pipeline's results to the 2_concatenate_and_filter\n pipeline's resources directory.\n"
	printf "\n The optional inflation and clustertype arguments select which MCL\n result directory to use; they do not affect the output dataset name.\n"
	printf "\n Usage: bash transition_1to2.sh <dataset> [inflation] [clustertype]\n"
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

if [ "$#" -gt 3 ]; then
	echo "[err] This scripts accepts maximum 3 arguments."
	usage
	exit 2
fi

DATASET=$1
if [ ! -d 1_cluster/results/${DATASET} ]; then
	echo "[err] Directory: 1_cluster/results/${DATASET} doesn't exist, first argument is wrong."
	usage
	exit 3
fi

INFLATION=${2:-1.8}
if [ ! -d 1_cluster/results/${DATASET}/mcl/$INFLATION ]; then
	echo "[err] Directory: 1_cluster/results/${DATASET}/mcl/$INFLATION doesn't exist, second argument is wrong or 1_cluster pipeline is net yet finished."
	usage
	exit 4
fi

CLUSTERTYPE=${3:-Normal}
if [ ! -d 1_cluster/results/${DATASET}/mcl/$INFLATION/$CLUSTERTYPE ]; then
	echo "[err] Directory: 1_cluster/results/${DATASET}/mcl/$INFLATION/$CLUSTERTYPE doesn't exist, third argument is wrong or 1_cluster pipeline is not yet finished."
	usage
	exit 5
fi

echo "[log] Creating directory $DATASET for dataset in 2_concatenate_and_filter pipeline ..."
mkdir -p 2_concatenate_and_filter/resources/${DATASET}/{mcl,eggnog}

echo "[log] Changing directory to 2_concatenate_and_filter/resources/${DATASET}/mcl ..."
cd 2_concatenate_and_filter/resources/${DATASET}/mcl

echo "[log] Symlinking sequences under 1_cluster/results/$DATASET/mcl/$INFLATION/$CLUSTERTYPE to 2_concatenate_and_filter pipeline's resources directory ..."
find ../../../../1_cluster/results/$DATASET/mcl/$INFLATION/$CLUSTERTYPE/ -iname "*.faa" -exec ln -s {} \;

echo "[log] Changing directory to 2_concatenate_and_filter/resources/${DATASET}/eggnog ..."
cd ../eggnog

echo "[log] Symlinking sequences under 1_cluster/results/$DATASET/cog_clusters/$CLUSTERTYPE to 2_concatenate_and_filter pipeline's resources directory ..."
find ../../../../1_cluster/results/$DATASET/cog_clusters/$CLUSTERTYPE/ -iname "*.faa" -exec ln -s {} \;

echo "[action required] Please update 2_concatenate_and_filter/config/config.yaml with the dataset name: $DATASET to perform the second step of the pipeline on it."

echo "[log] Script ended $(date)."

exit 0
