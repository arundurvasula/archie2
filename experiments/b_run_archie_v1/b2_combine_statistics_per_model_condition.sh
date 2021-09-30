#!/bin/bash

set -e

source ../export_vars.sh

SIMID=${SGE_TASK_ID}
if [ -z $SIMID ]; then
    SIMID=$1
fi

PARAMS=( $(head -n${SIMID} ../ms_testing_parameters.txt | tail -n1) )
MU=${PARAMS[0]}
R=${PARAMS[1]}
TESTDIR="$DATADIR/testing-data/msmodified.set-${SIMID}.mu-${MU}.r-${R}"

INPUTS=()
OUTPUT="$TESTDIR/archie_precision_recall_by_threshld.csv"
NREPS=50

if [ -e $OUTPUT ]; then
    echo "$OUTPUT already exists!"
    exit
fi

for REPL in `seq -f "%04g" 1 ${NREPS}`; do
    INPUT="$TESTDIR/$REPL/archie_statistics_per_snp.txt"
    if [ -e $INPUT ]; then
        INPUTS=( ${INPUTS[@]} $INPUT )
    fi
done


python3 combine_statistics.py -i ${INPUTS[@]} -o $OUTPUT 

