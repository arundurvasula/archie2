# This script simulates data to evaluate ARCHIE v1 on different model conditions, varying the mutation rate and recombination rate parameters. Except for these parameters, the simulations are the same as in the first ARCHIE paper (Durvasula and Sankararaman, 2019; also see https://github.com/sriramlab/ArchIE/tree/master/simulations). 

set -e

source ../export_vars.sh

SIMID=${SGE_TASK_ID}
if [ -z $SIMID ]; then
    SIMID=$1
fi

PARAMS=( $(head -n${SIMID} ../ms_testing_parameters.txt | tail -n1) )
MU=${PARAMS[0]}
R=${PARAMS[1]}

NTAR=100
NREF=100
NSITE=1000000
NREPS=100

MSLOG="ms.log"

TESTDIR="$DATADIR/testing-data/msmodified.set-${SIMID}.mu-${MU}.r-${R}"
mkdir -p $TESTDIR

for REPL in `seq -f "%04g" 1 ${NREPS}`; do
    echo "Processing replicate $REPL..."

    DIR="$TESTDIR/$REPL"
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
        if [ $REPL == "0001" ]; then
            SEED1=$(echo "3579 + $SIMID" | bc)
            SEED2=$(echo "27011 + $SIMID" | bc)
            SEED3=$(echo "59243 + $SIMID" | bc)
            echo "$SEED1 $SEED2 $SEED3" > seedms
        else
            cp $LASTDIR/seedms .
        fi

        echo "# Simulation model = $SIMID" > $MSLOG
        echo "# Recombination rate = $R" >> $MSLOG
        echo "# Mutation rate = $MU" >> $MSLOG
        echo "# Replicate = $REPL" >> $MSLOG

        SEEDS=( $(cat seedms) )
        echo "# seedms file before ms = ${SEEDS[@]}" >> $MSLOG

        bash $RUNSIMDIR/ms_model_durvasula2019.sh $NTAR $NREF $NSITE $MU $R >> $MSLOG

        SEEDS=( $(cat seedms) )
        echo "# seedms file after ms = ${SEEDS[@]}" >> $MSLOG

	echo "# done" >> $MSLOG
        echo "  ...done"
    fi

    LASTDIR="$DIR"

    echo "...done"
done

