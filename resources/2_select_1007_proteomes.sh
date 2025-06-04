#!/bin/bash

mkdir proteomes1007

cd proteomes1007

for f in `cat ../list_of_genomes_accession.txt`; do
	if [ -f ../protein_faa_reps/bacteria/${f}_protein.faa ]; then
		cp ../protein_faa_reps/bacteria/${f}_protein.faa ./
	else
		if [ -f ../protein_faa_reps/archaea/${f}_protein.faa ]; then
			cp ../protein_faa_reps/archaea/${f}_protein.faa ./
		else
			i=`echo $f | sed 's/GB_//' | sed 's/RS_//'`
			echo "$f not found, so downloading from NCBI using ID: $i"
			wget -q -O "${i}.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${i}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
			unzip -qq ${i}.zip
			if [ -f ncbi_dataset/data/${i}/protein.faa ]; then
				cp ncbi_dataset/data/${i}/protein.faa ./${i}_protein.faa
			else
				echo "no aa sequence available for ${i}, copying the na sequence"
				cp ncbi_dataset/data/${i}/${i}*_genomic.fna ./${i}_genomic.fna
			fi
			rm -rf ncbi_dataset
			rm README.md
			rm md5sum.txt
		fi
	fi
done

cd ..


echo "GCA_013046825.1 was published 2019 and is a draft genome of Turicibacter sanguinis, strain: DSM 14220"
echo "Maki et al 2020 sequences the same strain (albeit with different name: MOL361) and published a genome assembly https://pmc.ncbi.nlm.nih.gov/articles/PMC7303409/"
echo "with the accession ID: GCA_013046825.1"
echo "Let's use this insead of the missing genome GCA_004338625.1"

cd proteomes1007
i="GCA_013046825.1"
wget -q -O "${i}.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${i}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
unzip -qq ${i}.zip
if [ -f ncbi_dataset/data/${i}/protein.faa ]; then
	cp ncbi_dataset/data/${i}/protein.faa ./${i}_protein.faa
else
	echo "no aa sequence available for ${i}"
fi
rm -rf ncbi_dataset
rm README.md
rm md5sum.txt
