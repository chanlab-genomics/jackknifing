#!/bin/bash

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
ARRAY_TARGET=($DATA_DIR/Genomes_for_AFphylogeny/*)
# Specify the output directory
OUTPUT_DIR=$DATA_DIR/Genomes_for_AFphylogeny_red_40_26

# Create the output directory if it hasn't already been created
if [ ! -d $OUTPUT_DIR ] 
then
    mkdir $OUTPUT_DIR
fi


NCPUS=4
let CPUS_PER_TASK=2
let GROUP=$NCPUS/$CPUS_PER_TASK
let END=${#ARRAY_TARGET[@]}


for INDEX_OUTER in `seq 1 $NCPUS $END`; 
do
    let INNER_END=$INDEX_OUTER+$NCPUS-1
    for INDEX_INNER in `seq $INDEX_OUTER 1 $INNER_END`; 
    do
        if [ $INDEX_INNER -le $END ]
        then
            # This is where the sample command would go
            echo ${ARRAY_TARGET[$INDEX_INNER]}
        fi
    done

    # Wait after each grouping
    echo "wait"
done

time python3 ~/chanlab-genomics/jackknifing/jackknife.py --input_path /scratch/$PROJECT/$USER/Yeast/Genomes_for_AFphylogeny/AEG.fna --output_path /scratch/$PROJECT/$USER/Yeast/Genomes_for_AFphylogeny_red_40_26 --portion=40 --chunk_size=100 --threads=2