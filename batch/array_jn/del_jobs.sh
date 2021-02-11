#!/bin/bash

for JOB_ID in `seq 179506 1 179515`; 
do
    qdel $JOB_ID
done