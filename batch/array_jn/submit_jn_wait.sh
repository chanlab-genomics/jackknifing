#!/bin/bash

for SAMPLE_INDEX in `seq 26 1 27`;
do
    if [[ "$HOSTNAME" == *"gadi"* ]]; then
        qsub -N jn_ary_$i -P d85 -v -q normal -l wd SAMPLE_INDEX=$SAMPLE_INDEX ~/chanlab-genomics/jackknifing/batch/array_jn/jn_array.sh -o ~/chanlab-genomics/jackknifing/batch/jn_jobs/batch_out/jn_ary_$i.txt
    else
        qsub -N jn_ary_$i -A NCMAS-d85 SAMPLE_INDEX=$SAMPLE_INDEX ~/chanlab-genomics/jackknifing/batch/array_jn/jn_array.sh -o ~/chanlab-genomics/jackknifing/batch/jn_jobs/batch_out/jn_ary_$i.txt
    fi
done