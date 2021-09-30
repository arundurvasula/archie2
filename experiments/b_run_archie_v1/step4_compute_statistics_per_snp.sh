#!/bin/bash

set -e

source ../export_vars.sh

DIR=$1

INPUT="$DIR/archie_predictions_per_snp.txt"
OUTPUT="$DIR/archie_statistics_per_snp.txt"

if [ ! -e $INPUT ]; then
    echo "$INPUT does not exist!"
    exit 1
fi

DOAP=1
if [ -e $OUTPUT ]; then
    NLO=$(wc -l $OUTPUT | awk '{print $1}')
    if [ $NLO -eq 19 ]; then
        DOAP=0
    fi
fi

if [ $DOAP -eq 1 ]; then
    echo "  Computing ArchIE statistics per snp..."
    python3 compute_statistics_per_snp.py \
        -i $INPUT \
        -o $OUTPUT
    echo "  ...done"
fi

