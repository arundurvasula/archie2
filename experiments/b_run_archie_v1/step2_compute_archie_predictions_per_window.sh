#!/bin/bash

set -e

source ../export_vars.sh

DIR=$1

INPUT="$DIR/archie_features.txt"
OUTPUT="$DIR/archie_predictions.txt"

if [ ! -e $INPUT ]; then
    echo "$INPUT does not exist!"
    exit 1
fi

DOAP=1
if [ -e $OUTPUT ]; then
    NLI=$(wc -l $INPUT | awk '{print $1}')
    NLO=$(wc -l $OUTPUT | awk '{print $1}')
    if [ $[NLO-1] -eq $NLI ]; then
        DOAP=0
    fi
fi

if [ $DOAP -eq 1 ]; then
    echo "  Computing ArchIE predictions..."
    Rscript compute_archie_predictions.R $INPUT $OUTPUT
    echo "  ...done"
fi

