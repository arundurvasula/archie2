#!/bin/bash

for X in `seq 1 67`; do
qsub -cwd \
     -V \
     -N "a_sim_test_data" \
     -l h_data=32G,time=18:00:00,highp \
     -m e \
     -t $X \
     -M ekmolloy \
     -b y "bash a1_run_ms_for_testing.sh"
done

