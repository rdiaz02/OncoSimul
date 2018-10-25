#!/bin/bash

NRUNS=$1
INIT=1

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
    echo "Pass the number of processes as first argument. Optionally, starting as second" 
fi

if [[ $# -eq 2 ]]; then
    INIT=$2
    NRUNS=$(( $1 + $2 - 1 ))
fi

   
for (( run = $INIT; run <=${NRUNS}; run++ ))
do
    nohup R --vanilla --slave < test-manual-continuously.R &> test-MANUAL-continuously-$run.Rout &
done


## and search for "Failed" in tests
