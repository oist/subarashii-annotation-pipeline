#!/bin/bash

if [ $# -lt 3 ]; then
	echo "Too few arguments supplied."
	echo "Usage: $0 accession_ids.txt base_proteomes_dir output_dir"
	exit(1)
fi

IDS=$1
BASEDIR=$2
OUTDIR=$3

mkdir $OUTDIR

for f in `cat ${IDS} | awk -F "" '{print $4$5$6"/"$8$9$10"/"$11$12$13"/"$14$15$16"/"$4$5$6"_"$8$9$10$11$12$13$14$15$16"."$18"_genomic.fna.gz"}'`; do
	BASE=`basename -s _genomic.fna.gz $f`
	if [ -f ${BASEDIR}/$f ]; then
		cp ${BASEDIR}/$f ${OUTDIR}/
	else
		echo "$f not found, so downloading from NCBI using ID: $BASE"
		wget -q -O "${BASE}.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${BASE}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
		unzip -qq ${BASE}.zip
		if [ -f ncbi_dataset/data/${BASE}/${BASE}*genomic.fna ]; then
			cp ncbi_dataset/data/${BASE}/${BASE}*_genomic.fna ${OUTDIR}/${BASE}_genomic.fna
		else
			echo "no sequence found for ${BASE}"
		fi
		rm -rf ncbi_dataset
		rm README.md
		rm md5sum.txt
	fi
done
