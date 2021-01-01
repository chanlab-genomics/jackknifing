#!/bin/bash

for i in `seq 3 1 25`; 
do
    qsub -N jn_ary_$i -v INPUT_DIR=$(echo Genomes_for_AFphylogeny_red_40_$i) /home/s4430291/chanlab-genomics/jackknifing/jf_scripts/run_jellyfish.sh
done