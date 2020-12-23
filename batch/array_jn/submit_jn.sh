#!/bin/bash

for i in `seq 1 1 25`; 
do
    qsub -N jn_ary_$i -v OUTPUT_DIR="/30days/s4430291/Genomes_for_AFphylogeny_red_40_$i" /home/s4430291/chanlab-genomics/jackknifing/batch/array_jn/jn_array.sh
done