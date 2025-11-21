#!/bin/bash

usage() { 
	printf "\n$0 is a wrapper tool for EggNog mapper to run on a large number of files.\n
	Usage: $0 -i INPUT_DIR -e EGGNOG_DB_DIR -o OUTPUT_DIR [-c CPUs]\n
	INPUT_DIR: path to the directory containing the fasta files to annotate\n
	EGGNOG_DB_DIR: path to the directroy to the EggNog databases downloaded\n
	OUTPUT_DIR: path to the directory where the logfile and the annotations will be saved\n
	CPU: a number between 1 and 60\n" 1>&2; exit 1; 
}

CPU=1
while getopts "i:e:o:c:" o; do
	case "${o}" in
		i)
	    		INPUTDIR=${OPTARG}
	    		;;
		e)
	    		EGGNOGDB=${OPTARG}
	    		;;
		o)
	    		OUTPUT=${OPTARG}
	    		;;		
		c)
	    		CPU=${OPTARG}
	    		;;
		*)
	    		usage
	    		;;
    	esac
done
shift $((OPTIND-1))

if [ -z $INPUTDIR ] || [ -z $EGGNOGDB ] || [ -z $OUTPUT ]; then
	echo "Mandatory parameters are not set."
	usage
fi

if [ "$CPU" -lt "0" ] || [ "$CPU" -gt "60" ]; then
	CPU=1
fi

if [ ! -d ${INPUTDIR} ]; then
	echo "Directory $INPUTDIR does not exist."
	usage
fi

if [ ! -d $EGGNOGDB ]; then
	echo "EggNog DB is not present in $EGGNOGDB."
	usage
fi

mkdir -p ${OUTPUT}

MODE="diamond"

NOW=`date +'%Y-%m-%d_%H-%M'`
LOGFILE=`echo "$INPUTDIR-${NOW}.log" | sed 's/\//_/g' | sed 's/^\.//'`
echo $LOGFILE
echo "Per-genome EggNog annotation started on $INPUTDIR at $(date)" | tee ${OUTPUT}/${LOGFILE}
FILELIST=(`ls ${INPUTDIR}/`)

for f in ${FILELIST[@]}; do
	BASE=`basename $f`
	if [ ! -f ${OUTPUT}/${BASE}.emapper.annotations ]; then
		echo "Starting $BASE ..." | tee -a ${OUTPUT}/${LOGFILE}
		emapper.py -i ${INPUTDIR}/${f} --cpu ${CPU} --data_dir ${EGGNOGDB} -m ${MODE} --output ${OUTPUT}/${BASE}.emapper --override 2>&1 | tee -a ${OUTPUT}/${LOGFILE}
		echo "Done." | tee -a ${OUTPUT}/${LOGFILE}
	else
		echo "Skipping $BASE (already annotated)" | tee -a ${OUTPUT}/${LOGFILE}
	fi
done
echo "Finished at $(date)" | tee -a ${OUTPUT}/${LOGFILE}
