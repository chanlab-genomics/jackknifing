#!/bin/bash

# tar -czf /scratch/$PROJECT/$USER/Yeast/Genomes_for_AFphylogeny_red_40_26_D2S.tz.gz /scratch/$PROJECT/$USER/Yeast/Genomes_for_AFphylogeny_red_40_26_D2S

for SAMPLE_INDEX in `seq 27 1 50`;
do
    echo "Compressing sample ${SAMPLE_INDEX}"
    # Use the -P option to use absolute path names
    tar -czf /scratch/$PROJECT/$USER/Yeast/D2S_archive/Genomes_for_AFphylogeny_red_40_${SAMPLE_INDEX}_D2S.tz.gz /scratch/$PROJECT/$USER/Yeast/Genomes_for_AFphylogeny_red_40_${SAMPLE_INDEX}_D2S
done