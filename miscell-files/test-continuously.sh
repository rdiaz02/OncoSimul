#!/bin/bash

NRUNS=$1

if [[ $# -ne 1 ]]; then
   echo "Pass the number of processes as argument"
fi

   
   
for (( run = 1; run <=${NRUNS}; run++ ))
do
    nohup R --vanilla --slave < test-continuously.R &> test-continuously-$run.Rout &
done

