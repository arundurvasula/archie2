#!/bin/bash

for X in `seq 1 49`; do
qsub -cwd \
     -V \
     -N "b_run_archie_v1" \
     -l h_data=32G,time=8:00:00,highp \
     -m e \
     -t $X \
     -M ekmolloy \
     -b y "bash b1_run_archie_on_ms_testing.sh"
done

