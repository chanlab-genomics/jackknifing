#!/bin/bash

for SAMPLE_INDEX in `seq 56 1 56`;
do
    if [[ "$HOSTNAME" == *"gadi"* ]]; then
        qsub -N jn_ary_$SAMPLE_INDEX -P d85 -q normal -l wd -v SAMPLE_INDEX=$SAMPLE_INDEX -o ~/chanlab-genomics/jackknifing/batch/jn_jobs/batch_out/jn_ary_$SAMPLE_INDEX.txt ~/chanlab-genomics/jackknifing/batch/array_jn/jn_array.sh
    else
        qsub -N jn_ary_$SAMPLE_INDEX -A NCMAS-d85 -v SAMPLE_INDEX=$SAMPLE_INDEX -o ~/chanlab-genomics/jackknifing/batch/jn_jobs/batch_out/jn_ary_$SAMPLE_INDEX.txt ~/chanlab-genomics/jackknifing/batch/array_jn/jn_array.sh
    fi

    sleep 1
done