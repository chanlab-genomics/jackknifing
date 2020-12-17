#!/bin/bash
#
#PBS -N jn_py
#-#PBS -l nodes=1
#PBS -l select=1:ncpus=8:mem=50GB
#PBS -l walltime=00:30:00

#PBS -j oe
#PBS -o /gpfs1/homes/s4430291/chanlab-genomics/jackknifing/example_out/jn_py_bio.txt

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

module load python
time python3 /gpfs1/homes/s4430291/chanlab-genomics/jackknifing/jackknife.py --input_path /30days/s4430291/Slin_CCMP2456_py/S.linucheae_CCMP2456.genome.fasta --output_path /30days/s4430291/Slin_CCMP2456_py_out --portion=40 --chunk_size=100 --threads=8

echo "time finished "$DATE

#Refer to the table of job environment variables section of the RCC PBS Pro User Guide