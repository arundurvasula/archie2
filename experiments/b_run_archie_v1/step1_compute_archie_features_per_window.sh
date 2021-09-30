#!/bin/bash

set -e

source ../export_vars.sh

DIR=$1
WINDOW_SIZE=$2 # 50000
STEP_SIZE=$3   # 10000

AFOUT="$DIR/archie_features.txt"
AFLOG="$DIR/archie_features.log"

NHAP=$(wc -l $DIR/out.ADMIXED.ind | awk '{print $1}')
CHR=$(head -n1 $DIR/out.snp | sed 's/:/ /g' | awk '{print $1}')
START=$(head -n1 $DIR/out.snp | sed 's/:/ /g' | awk '{print $2}')
END=$(tail -n1 $DIR/out.snp | sed 's/:/ /g' | awk '{print $2}')

# Check ArchIE features completed
DOAF=1
if [ -e $AFOUT ]; then
    WINDOWS=( $(seq $START $STEP_SIZE $END) )
    NLINES=$( echo "${#WINDOWS[@]} * $NHAP" | bc )
    if [ $NLINES -eq $(wc -l $AFOUT | awk '{print $1}') ]; then
        DOAF=0
    fi
fi

if [ $DOAF -eq 1 ]; then
    echo "  Computing ArchIE features..."
    if [ $(wc -l $DIR/out.snp | awk '{print $1}') -gt 0 ]; then
        python3 $COMPUTE_AF \
            -s $DIR/out.snp \
            -i $DIR/out.ADMIXED.ind \
            -a $DIR/out.ADMIXED.geno \
            -r $DIR/out.1.geno \
            -c $CHR \
            -b $START \
            -e $END \
            -w $WINDOW_SIZE \
            -z $STEP_SIZE 1> $AFOUT 2> $AFLOG
    fi
    echo "  ...done"
fi

