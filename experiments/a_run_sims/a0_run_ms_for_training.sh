# This script simulates training data as described in the first ARCHIE paper (Durvasula and Sankararaman, 2019). It is not necessary to run this script to perform the experiments that compare ARCHIE v1 on model conditions with different mutation rate and recombination rate parameters.

set -e

source ../export_vars.sh

SIMID=${SGE_TASK_ID}
if [ -z $SIMID ]; then
   SIMID=$1
fi

PARAMS=( $(head -n${SIMID} ../ms_training_parameters.txt | tail -n1) )
MU=${PARAMS[0]}
R=${PARAMS[1]}

NTAR=100
NREF=100
NSITE=50000
NREPS=10000

MSLOG="ms.log"

TRAINDIR="$DATADIR/training-data/msmodified.set-${SIMID}.mu-${MU}.r-${R}"
mkdir -p $TRAINDIR

for REPL in `seq -f "%05g" 1 ${NREPS}`; do
    echo "Processing replicate $REPL..."

    DIR="$TRAINDIR/$REPL"
    mkdir -p $DIR
    cd $DIR

    DOMS=1
    if [ -e $MSLOG ]; then
         if [ ! -z $(grep "Admixture proportion" $MSLOG | awk '{print $2}') ]; then
             DOMS=0
         fi
    fi

    if [ $DOMS -eq 1 ]; then
        echo "  Simulating data with ms..."
        if [ $REPL == "00001" ]; then
            echo "3579 27011 59243" > seedms
        else
            cp $LASTDIR/seedms .
        fi

        echo "# Simulation model = $SIMID" > $MSLOG
        echo "# Recombination rate = $R" >> $MSLOG
        echo "# Mutation rate = $MU" >> $MSLOG
        echo "# Replicate = $REPL" >> $MSLOG

        SEEDS=( $(cat seedms) )
        echo "# seedms file before = ${SEEDS[@]}" >> $MSLOG

        bash $RUNSIMDIR/ms_model_durvasula2019.sh $NTAR $NREF $NSITE $MU $R >> $MSLOG

        SEEDS=( $(cat seedms) )
        echo "# seedms file after = ${SEEDS[@]}" >> $MSLOG

        echo "# done" >> $MSLOG
        echo "  ...done"
    fi

    LASTDIR="$DIR"

    echo "...done"
done

