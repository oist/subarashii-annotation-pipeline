#!/bin/bash

usage(){
	echo "This program will summarize the resource usage of EggNog annotation analysis initiated by a snakemake job on a SLURM cluster"
	echo
	echo "Usage: $0 <path to snakemake logfile> <output file>"
}

if [ $# -lt 2 ]; then
	echo "Error: not enough paramters"
	usage
	exit 1
fi
if [ $# -gt 2 ]; then
	echo "Error: too many parameters"
	usage
	exit 2
fi

SNAKEMAKE_LOGFILE=$1
OUTFILE=$2

printf "jobid,dataset,filename,numseqs,seqlenavg,seqlenmedian,seqlenmin,seqlenmax,above1k,above10k,maxmem,cpus,time,hits\n" > $OUTFILE

COUNTER=0
SUM_USED_MEM=0
for id in `grep "SLURM jobid" $SNAKEMAKE_LOGFILE | cut -d " " -f 9`; do
	#echo "ID: $id"
	INFO=`sacct -j $id --format="JobId%30,JobName,NCPUS,ElapsedRaw,StdOut%200,MaxRSS" | tail -n +3 | head -n 1 | sed -E 's/\s+/ /g' | sed -E 's/^\s+//' | sed -E 's/\s+$//'`
	MAX_MEM=`sacct -j $id --format="MaxRSS" | tail -n 1 | sed -E 's/\s+/ /g' | sed -E 's/^\s+//' | sed -E 's/\s+$//' | sed 's/K//'`
	echo "Info: $INFO ${MAX_MEM}K"
	CPUS=`echo $INFO | cut -d " " -f 3`
	#echo "CPUs: $CPUS"
	ELAPSED=`echo $INFO | cut -d " " -f 4`
	#echo "Elapsed: $ELAPSED"
	LOGFILE=`echo $INFO | cut -d " " -f 5`
	#echo "Logfile: $LOGFILE"
	if [ -z "${MAX_MEM}" ]; then
		MAX_MEM=0
	fi
	#echo "Max mem: $MAX_MEM"
	DATASET=`echo $LOGFILE | sed -E 's/^.+rule_eggnog_mapper\/([a-zA-Z_0-9]+)_([0-9]+)\/\%j.log$/\1 \2/' | cut -d " " -f 1`
	#echo "Dataset: $DATASET"
	PART=`echo $LOGFILE | sed -E 's/^.+rule_eggnog_mapper\/([a-zA-Z_0-9]+)_([0-9]+)\/\%j.log$/\1 \2/' | cut -d " " -f 2`
	#echo "PART: $PART"
	# max mem is only set for finished jobs, so this is a filter for already finished jobsq
	if [ $MAX_MEM -gt 0 ]; then
		SUM_USED_MEM=$((SUM_USED_MEM+MAX_MEM))
		COUNTER=$((COUNTER+1))

		# calculate sequence length
		python3 ./seq_length.py -i ../resources/${DATASET}/all_genomes_protein.part${PART}.faa -t fasta -o part${PART}.seqlen
		NUMSEQS=`wc -l part${PART}.seqlen | sed -E 's/\s+/ /g' | cut -d " " -f 1`
		SEQLEN=`sort -k2 part${PART}.seqlen | awk 'BEGIN{sum=0; count=0; min=1e9; max=-1e9; above1k=0; above10k=0} {lengths[NR]=$2; if ($2 < min) { min = $2 }; if ($2 > max) { max = $2 }; sum+=$2; count+=1; if ( $2 > 1000) {above1k+=1}; if ($2 > 10000) {above10k+=1} } END{if (NR % 2) { print sum/NR, lengths[(NR+1) / 2],  min, max, above1k, above10k } else {print sum/NR, ((lengths[(NR/2)] + lengths[(NR/2)+1]) / 2.0), min, max, above1k, above10k } }'`
		
		AVGSEQLEN=`echo $SEQLEN | cut -d " " -f 1`
		MEDIANSEQLEN=`echo $SEQLEN | cut -d " " -f 2`
		MINSEQLEN=`echo $SEQLEN | cut -d " " -f 3`
		MAXSEQLEN=`echo $SEQLEN | cut -d " " -f 4`
		ABOVE1K=`echo $SEQLEN | cut -d " " -f 5`
		ABOVE10K=`echo $SEQLEN | cut -d " " -f 6`

		ELAPSED_HOURS=`echo $ELAPSED | awk '{print $1 / 60 /60}'`

		HITS=`wc -l ../results/${DATASET}/eggnog/all_genomes_protein.part${PART}.emapper.annotations | sed -E 's/\s+/ /g' | cut -d " " -f 1`
		MAX_MEM_GB=`echo $MAX_MEM | awk '{print $1 / 1024 / 1024}'`
		echo "Max mem: $MAX_MEM_GB"
		printf "$id,${DATASET},all_genomes_protein.part${PART}.faa,$NUMSEQS,$AVGSEQLEN,$MEDIANSEQLEN,$MINSEQLEN,$MAXSEQLEN,${ABOVE1K},${ABOVE10K},$MAX_MEM_GB,$CPUS,$ELAPSED_HOURS,$HITS\n" >> $OUTFILE
	fi
done

#AVG_MEM_USE=awk 'BEGIN{sum=0; count=0} {sum+=$1; count+=1} END{print sum/count /1024 /1024 "GB"}'



