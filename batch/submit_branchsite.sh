#!/bin/bash

ortholist=multi6
walltime=6:00:00

for i in {1..500}
do
    for model in Alt7 # Models: Alt1 Alt2 Alt3 Alt4 Alt5 Alt6 Alt7 Alt8 Alt9 Null1 Null2 Null3
    do
       orthoname=$(awk 'NR=='${i}'' orthogroup_lists/branch_site_orthogroups_${ortholist}.txt | awk '{print $NF}' FS=/ )
       jobname=$(echo ${orthoname}_${model}_${ortholist})
       qsub -l walltime=$walltime -N $jobname -v ORTHOLIST=$ortholist,MODEL=$model,ORTHONAME=$orthoname,JOBNAME=$jobname branchsite_all.pbs
    done
done
