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

NREPS=100

for REPL in `seq -f "%04g" 1 ${NREPS}`; do
    echo "Processing replicate $REPL..."

    DIR="$TESTDIR/$REPL"

    # Compute ArchIE features
    bash step1_compute_archie_features_per_window.sh $DIR 50000 10000

    # Compute ArchIE predictions in each window using trained model
    bash step2_compute_archie_predictions_per_window.sh $DIR

    # Average prediction per SNP
    bash step3_average_predictions_per_snp.sh $DIR 

    # Compute TP, FP, TN, FN
    bash step4_compute_statistics_per_snp.sh $DIR

    echo "...done"
done

