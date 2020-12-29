#!/bin/bash

NRUNS=$1
INIT=1

## Actually, when left to run unattended, multiple, tmux, on detached,
## all except one seem to fail in the vignette part.
## Probably because they are writing to the same files.

if [ $# -ne 1 ] && [ $# -ne 2 ]; then
    echo "Pass the number of processes as first argument. Optionally, starting as second" 
fi

if [[ $# -eq 2 ]]; then
    INIT=$2
    NRUNS=$(( $1 + $2 - 1 ))
fi

   
for (( run = $INIT; run <=${NRUNS}; run++ ))
do
    nohup R --vanilla --slave < test-continuously.R &> test-continuously-$run.Rout &
done


## grep for FAILURE, "Test failures", FAIL, Failed
