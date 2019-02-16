#!/bin/bash


NRUNS=$1
SCRIPTFULL=$2
SCRIPT="${SCRIPTFULL%.*}"

if [[ $# -ne 2 ]]; then
   echo "Pass the number of processes and script name as argument"
fi

randomstring () {
    cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w ${1:-32} | head -n 1
}

THISEXEC=$(randomstring)

for (( run = 1; run <=${NRUNS}; run++ ))
do
    nohup R --vanilla --slave < $SCRIPTFULL &> $SCRIPT-$THISEXEC-$run.Rout &
done

