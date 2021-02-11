#!/bin/bash

for SAMPLE_INDEX in `seq 101 1 115`;
do

    if [[ "$HOSTNAME" == *"gadi"* ]]; then
        qsub -N D2S_jf_$SAMPLE_INDEX -P d85 -q normal -l wd -v INPUT_DIR=$(echo Genomes_for_AFphylogeny_red_40_$SAMPLE_INDEX) -o ~/chanlab-genomics/jackknifing/batch/jellyfish_jobs/batch_out/jf_ary_$SAMPLE_INDEX.txt ~/chanlab-genomics/jackknifing/batch/jellyfish_jobs/run_jellyfish.sh
    else
        qsub -N D2S_jf_$SAMPLE_INDEX -A NCMAS-d85 INPUT_DIR=$(echo Genomes_for_AFphylogeny_red_40_$SAMPLE_INDEX) -o ~/chanlab-genomics/jackknifing/batch/jellyfish_jobs/batch_out/jf_ary_$SAMPLE_INDEX.txt ~/chanlab-genomics/jackknifing/batch/jellyfish_jobs/run_jellyfish.sh
    fi

    sleep 1
done