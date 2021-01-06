#!/bin/bash

# TARGET_DIR=$(echo /30days/$USER/Genomes_for_AFphylogeny)
TARGET_DIR=$(echo /scratch/$PROJECT/$USER/Yeast/Genomes_for_AFphylogeny)

OUTPUT_DIR_PREFIX=$(echo /scratch/$PROJECT/$USER/Yeast)

for i in `seq 1 1 1`; 
do
    qsub -N jn_ary_$i -v OUTPUT_DIR=$(echo $OUTPUT_DIR_PREFIX/Genomes_for_AFphylogeny_red_40_$i) ~/chanlab-genomics/jackknifing/batch/array_jn/jn_array_gadi.sh
done