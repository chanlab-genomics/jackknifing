#!/bin/bash

for JOB_ID in `seq 510889 1 510908`; 
do
    qdel $JOB_ID
done