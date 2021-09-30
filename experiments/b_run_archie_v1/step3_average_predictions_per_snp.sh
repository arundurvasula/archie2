#!/bin/bash

set -e

source ../export_vars.sh

DIR=$1

INPUT="$DIR/archie_predictions.txt"
SNPS="$DIR/out.snp"
ANCS="$DIR/out.ADMIXED.anc"
OUTPUT="$DIR/archie_predictions_per_snp.txt"

if [ ! -e $INPUT ]; then
    echo "$INPUT does not exist!"
    exit 1
fi

DOAP=1
if [ -e $OUTPUT ]; then
    NLI=$(wc -l $SNPS | awk '{print $1}')
    NLI=$(echo "($NLI * 100) + 1" | bc)
    NLO=$(wc -l $OUTPUT | awk '{print $1}')
    if [ $NLO -eq $NLI ]; then
        DOAP=0
    fi
fi

if [ $DOAP -eq 1 ]; then
    echo "  Averaging ArchIE predictions per snp..."
    python3 average_predictions_per_snp.py \
        -i $INPUT \
        -s $SNPS \
        -a $ANCS \
        -o $OUTPUT
    echo "  ...done"
fi

