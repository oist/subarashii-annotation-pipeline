#!/bin/bash

echo "This script is not needed, as we'll use the aa sequences."
#exit 0

mkdir genomes1007

cd genomes1007

for f in `cat ../list_of_genomes_accession.txt | awk -F "" '{print $4$5$6"/"$8$9$10"/"$11$12$13"/"$14$15$16"/"$4$5$6"_"$8$9$10$11$12$13$14$15$16"."$18"_genomic.fna.gz"}'`; do
	BASE=`basename -s _genomic.fna.gz $f`
	if [ -f ../gtdb_genomes_reps_r202/$f ]; then
		cp ../gtdb_genomes_reps_r202/$f ./
	else
		echo "$f not found, so downloading from NCBI using ID: $BASE"
		wget -q -O "${BASE}.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${BASE}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
		unzip -qq ${BASE}.zip
		if [ -f ncbi_dataset/data/${BASE}/${BASE}*genomic.fna ]; then
			cp ncbi_dataset/data/${BASE}/${BASE}*_genomic.fna ./${BASE}_genomic.fna
		else
			echo "no sequence found for ${BASE}"
		fi
		rm -rf ncbi_dataset
		rm README.md
		rm md5sum.txt
	fi
done

cd ..

echo "GCA_013046825.1 was published 2019 and is a draft genome of Turicibacter sanguinis, strain: DSM 14220"
echo "Maki et al 2020 sequences the same strain (albeit with different name: MOL361) and published a genome assembly https://pmc.ncbi.nlm.nih.gov/articles/PMC7303409/"
echo "with the accession ID: GCA_013046825.1"
echo "Let's use this insead of the missing genome GCA_004338625.1"

cd genomes1007
i="GCA_013046825.1"
wget -q -O "${i}.zip" "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${i}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
unzip -qq ${i}.zip
if [ -f ncbi_dataset/data/${i}/${i}*genomic.fna ]; then
        cp ncbi_dataset/data/${i}/${i}*genomic.fna ./${i}_genomic.fna
else
        echo "no sequence available for ${i}"
fi
rm -rf ncbi_dataset
rm README.md
rm md5sum.txt

