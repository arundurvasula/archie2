#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load gcc/4.9.3
module load python/3.7.2
module load R/3.4.0

export HOMEDIR="$HOME/project-sriram/archie2/experiments"
export EXTDIR="$HOMEDIR/external"
export RUNSIMDIR="$HOMEDIR/a_run_sims"
export RUNARCHIEDIR="$HOMEDIR/b_run_archie"
export DATADIR="/u/scratch/e/ekmolloy/archie2/experiments/data"

export COMPUTE_AF_FOR_TRAINING="$EXTDIR/ArchIE/simulations/calc_stats_ms.py"
export COMPUTE_AF="$EXTDIR/ArchIE/data/calc_stats_window_data.py"

export PATH="$PYBINDIR:$PATH"
export PATH="$EXTDIR/ArchIE/msmodified:$PATH"  # ms lives here
export PYTHONPATH="$HOME/.local/lib/python3.7/site-packages:$PYTHONPATH"

