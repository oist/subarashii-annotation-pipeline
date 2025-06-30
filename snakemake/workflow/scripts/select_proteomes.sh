#!/bin/bash

if [ $# -lt 3 ]; then
	echo "Too few arguments supplied."
	echo "Usage: $0 accession_ids.txt base_proteomes_dir output_dir"
	exit 1
fi

IDS=$1
BASEDIR=$2
OUTDIR=$3

mkdir ${OUTDIR}

for f in `cat ${IDS}`; do
	if [ -f ${BASEDIR}/bacteria/${f}_protein.faa ]; then
		cp ${BASEDIR}/bacteria/${f}_protein.faa ./
	else
		if [ -f ${BASEDIR}/archaea/${f}_protein.faa ]; then
			cp ${BASEDIR}/archaea/${f}_protein.faa ./
		else
			i=`echo $f | sed 's/GB_//' | sed 's/RS_//'`
			echo "$f not found, so downloading from NCBI using ID: $i"
			wget -q -O "${i}.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${i}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
			unzip -qq ${i}.zip
			if [ -f ncbi_dataset/data/${i}/protein.faa ]; then
				cp ncbi_dataset/data/${i}/protein.faa ${OUTDIR}/${i}_protein.faa
			else
				echo "no aa sequence available for ${i}, copying the na sequence"
				cp ncbi_dataset/data/${i}/${i}*_genomic.fna ${OUTDIR}/${i}_genomic.fna
			fi
			rm -rf ncbi_dataset
			rm README.md
			rm md5sum.txt
		fi
	fi
done
