#!/bin/bash
#
#-#PBS -l nodes=1
#PBS -l select=1:ncpus=6:mem=30GB
#PBS -l walltime=00:20:00

# Make this a job array
#PBS -J 0-196

#PBS -j oe

#CHANGE THIS TO YOUR UQ-FACULTY-SCHOOL group name. 
#USE the groups command to find out your exact group name. 
#PBS -A NCMAS-d85

#PBS -l select=1
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

# Make a array of folders of the target directory
ARRAY_TARGET=(/30days/s4430291/Genomes_for_AFphylogeny/*.fna)

# Specify the output directory
# OUTPUT_DIR=/30days/s4430291/Genomes_for_AFphylogeny_red_40_$1
echo OUTPUT_DIR: $OUTPUT_DIR

if [ ! -d $OUTPUT_DIR ] 
then
    mkdir $OUTPUT_DIR
fi

module load python
python3 /gpfs1/homes/s4430291/chanlab-genomics/jackknifing/jackknife.py --input_path ${ARRAY_TARGET[$PBS_ARRAY_INDEX]} --output_path $OUTPUT_DIR --portion=40 --chunk_size=100 --threads=$NCPUS

echo "time finished "$DATE