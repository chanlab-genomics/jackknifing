#!/bin/bash
#PBS -l mem=180GB,walltime=6:00:00,ncpus=12
#PBS -l jobfs=20GB
#PBS -k oe
#PBS -j oe

#PBS -A NCMAS-d85

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

set -o errexit
cd $PBS_O_WORKDIR

# Check which machine we are on and adjust the directory to the yeast genome 
# data accordingly
if [[ "$HOSTNAME" == *"gadi"* ]]; then
    DATA_DIR=/scratch/$PROJECT/$USER/Yeast
    module load python3/3.7.4
    module load python2/
    module load jellyfish/2.3.0
else
    DATA_DIR=/30days/$USER
    module load python
    module load jellyfish/2.2.10
fi

k=21

ARRAY_TARGET=($DATA_DIR/$INPUT_DIR/*.fna)

s=10000000000

# NOTE: I changed jellyfish dump -ct $file.$k.jf_0 to jellyfish dump -ct $file.$k.jf
# in response to the following error in pbs output
## File: /30days/s4430291/Genomes_for_AFphylogeny_red_40/AEG_40.fna
# Failed to open input file '/30days/s4430291/Genomes_for_AFphylogeny_red_40/AEG_40.fna.21.jf_0'


let END=${#ARRAY_TARGET[@]}-1
let CPUS_PER_TASK=4
let GROUP=$NCPUS/$CPUS_PER_TASK

for INDEX_OUTER in `seq 0 $GROUP $END`; 
do
    let INNER_END=$INDEX_OUTER+$GROUP-1
    for INDEX_INNER in `seq $INDEX_OUTER 1 $INNER_END`; 
    do
        if [ $INDEX_INNER -le $END ]
        then
            # This is where the sample command would go
            # echo Starting ${ARRAY_TARGET[$INDEX_INNER]}
            # python3 ~/chanlab-genomics/jackknifing/jackknife.py --input_path ${ARRAY_TARGET[$INDEX_INNER]} --output_path $OUTPUT_DIR --portion=40 --chunk_size=100 --threads=$CPUS_PER_TASK &
        
            # Run a given line in the array
            file=${ARRAY_TARGET[$INDEX_INNER]}
            jellyfish count -m $k -s $s -t $CPUS_PER_TASK -o $file.$k.jf $file && \
            jellyfish dump -ct $file.$k.jf | sort -k1,1 | python2 ~/chanlab-genomics/jackknifing/jf_scripts/Kmers_2_NumbericRepresentation.py -o $file.${k}mer.nkc.gz && \
            python2 ~/chanlab-genomics/jackknifing/jf_scripts/Composition_of_InputSeqs.py --fasta $file --freq $file.CharFreq && \
            touch $file.done &
        fi
    done

    # Wait after each grouping
    wait
done

DATE=$(date +"%d/%m/%Y %H:%M")
echo "time finished "$DATE