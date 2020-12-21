#!/bin/bash
#
#PBS -N test_array
#-#PBS -l nodes=1
#PBS -l select=1:ncpus=1:mem=1GB
#PBS -l walltime=00:00:10

# Make this a job array
#PBS -J 0-5

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

cd /home/s4430291/chanlab-genomics/jackknifing/example_out/tmp
pwd

echo arrary num: $PBS_ARRAY_INDEX

echo "time finished "$DATE