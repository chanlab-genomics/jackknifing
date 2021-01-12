#!/bin/bash

#PBS -l ncpus=8,mem=5GB
#PBS -l walltime=00:20:00

#PBS -j oe

DATE=$(date +"%d/%m/%Y %H:%M")
echo "time started  "$DATE
echo ------------------------------------------------------
echo -n 'Job is running on the following nodes '; cat $PBS_NODEFILE
echo ------------------------------------------------------
echo PBS: qsub is running on $PBS_O_HOST
echo PBS: originating queue is $PBS_O_QUEUE
echo PBS: executing queue is $PBS_QUEUE
echo PBS: working directory is $PBS_O_WORKDIR
echo PBS: execution mode is $PBS_ENVIRONMENT
echo PBS: job identifier is $PBS_JOBID
echo PBS: job name is $PBS_JOBNAME
echo PBS: node file is $PBS_NODEFILE
echo PBS: current home directory is $PBS_O_HOME
echo PBS: PATH = $PBS_O_PATH
echo ------------------------------------------------------
export TIMEFORMAT="%E sec"

cd $PBS_O_WORKDIR
pwd

# Check which machine we are on and adjust the directory to the yeast genome 
# data accordingly
if [[ "$HOSTNAME" == *"gadi"* ]]; then
    DATA_DIR=/scratch/$PROJECT/$USER/Yeast
    module load python3/3.7.4
else
    DATA_DIR=/30days/$USER
    module load python
fi

# Make a array of folders of the target directory
ARRAY_TARGET=($DATA_DIR/Genomes_for_AFphylogeny/*.fna)
# Specify the output directory
OUTPUT_DIR=$DATA_DIR/Genomes_for_AFphylogeny_red_40_$SAMPLE_INDEX

# Create the output directory if it hasn't already been created
if [ ! -d $OUTPUT_DIR ] 
then
    mkdir $OUTPUT_DIR
fi

let END=${#ARRAY_TARGET[@]}-1
let CPUS_PER_TASK=2
let GROUP=$NCPUS/$CPUS_PER_TASK

for INDEX_OUTER in `seq 0 $GROUP $END`; 
do
    let INNER_END=$INDEX_OUTER+$GROUP-1
    for INDEX_INNER in `seq $INDEX_OUTER 1 $INNER_END`; 
    do
        if [ $INDEX_INNER -le $END ]
        then
            # This is where the sample command would go
            echo Starting ${ARRAY_TARGET[$INDEX_INNER]}
            python3 ~/chanlab-genomics/jackknifing/jackknife.py --input_path ${ARRAY_TARGET[$INDEX_INNER]} --output_path $OUTPUT_DIR --portion=40 --chunk_size=100 --threads=$CPUS_PER_TASK &
        fi
    done

    # Wait after each grouping
    wait
done

DATE=$(date +"%d/%m/%Y %H:%M")
echo "time finished "$DATE