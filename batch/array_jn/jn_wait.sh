
# Make a array of folders of the target directory
ARRAY_TARGET=(/30days/s4430291/Genomes_for_AFphylogeny/*)

# Specify the output directory
OUTPUT_DIR=/30days/s4430291/Genomes_for_AFphylogeny_red_40_26

if [ ! -d $OUTPUT_DIR ] 
then
    mkdir $OUTPUT_DIR
fi

module load python

NCPUS=4
let END=${#ARRAY_TARGET[@]}

for INDEX_OUTER in `seq 1 $NCPUS $END`; 
do
    let INNER_END=$INDEX_OUTER+$NCPUS-1
    for INDEX_INNER in `seq $INDEX_OUTER 1 $INNER_END`; 
    do
        if [ $INDEX_INNER -le $END ]
        then
            # This is where the sample command would go
            echo ${ARRAY_TARGET[$INDEX_INNER]}
        fi
    done

    # Wait after each grouping
    echo "wait"
done