#!/bin/bash

if [[ "$HOSTNAME" == *"gadi"* ]]; then
    DATA_IN_DIR=/scratch/$PROJECT/$USER/Yeast
    DATA_OUT_DIR=/scratch/$PROJECT/$USER/Yeast
else
    DATA_OUT_DIR=/30days/$USER
    DATA_OUT_DIR=/90days/$USER
fi

let INDEX=0
# python3 calc_d2s/create_d2s_jobs.py --data_input_path $DATA_IN_DIR/Genomes_for_AFphylogeny_red_40_${INDEX} --data_output_path $DATA_OUT_DIR/Genomes_for_AFphylogeny_red_40_${INDEX}_D2S --temp F --submit F --dry_run F --index=${INDEX}
python3 calc_d2s/create_d2s_jobs.py --data_input_path /30days/s4430291/Genomes_for_AFphylogeny --data_output_path /90days/s4430291/Genomes_for_AFphylogeny --temp F --submit T --dry_run F --index=${INDEX}