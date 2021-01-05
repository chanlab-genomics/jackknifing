#!/bin/bash

let END=99
let GROUP=10

for INDEX_OUTER in `seq 1 $GROUP $END`; 
do
    let INNER_END=$INDEX_OUTER+$GROUP-1
    for INDEX_INNER in `seq $INDEX_OUTER 1 $INNER_END`; 
    do
        if [ $INDEX_INNER -le $END ]
        then
            # This is where the sample command would go
            echo $INDEX_INNER
        fi
    done

    # Wait after each grouping
    echo "wait"
done