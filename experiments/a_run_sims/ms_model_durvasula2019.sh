# This script simulates data based on Figure 1 in the first ARCHIE paper (Durvasula and Sankararaman, 2019), although note the bottlenecking for Neanderthal population. This script is based on https://github.com/sriramlab/ArchIE/tree/master/simulations 

set -e

# Set variable model parameters
NCEU=$1   # Sample size of target population (CEU)
          # Original value: 100
NYRI=$2   # Sample size of reference population (YRI)
          # Original value: 100
NSITE=$3  # No. of sites per locus (= window size in ArchIE) 
          # Original value: 50000
MU=$4     # Mutation rate (i.e. event per site per allele)
          # Original value: 0.0000000125
R=$5      # Recombination rate (i.e. event per site per allele)
          # Original value: 0.00000001

# Set fixed model parameters
NREPL=1       # Number of independent replicates (loci) to simulate
PSIZE=10000   # Effective population size (no. of alleles in each population)
NPOP=4        # No. of populations
NNEA=1        # Neanderthal sample size
NDEN=1        # Denisovan sample size
AFRAC="0.02"  # Admixture fraction going forward in time
BFRAC="0.01"  # Bottleneck fraction for Neanderthals
T_A=2000      # Time that admixture occurs
T_Y2C=2500    # Time that YRI joins CEU (Modern Humans)
T_NBS=6000    # Time that Neanderthal bottleneck start 
T_NBE=6120    # Time that Neanderthal bottleneck end
T_D2N=7000    # Time that Denisovan joins Neanderthal (Archaics)
T_0=12000     # Time that Archaics join Modern Humans

# Convert model parameters for use with ms
NSAMP=$(echo "$NCEU + $NYRI + $NNEA + $NDEN" | bc)
THETA=$(echo "4 * $PSIZE * $MU * $NSITE" | bc)  # Expected no. of new mutations in locus
RHO=$(echo "4 * $PSIZE * $R * $NSITE" | bc)     # Expected no. of recombination breakpoints in locus
AFRAC=$(echo "1 - $AFRAC" | bc)                 # Subtract from 1, because going backward in time
T_A=$(echo "$T_A / (4 * $PSIZE)" | bc -l)
T_Y2C=$(echo "$T_Y2C / (4 * $PSIZE)" | bc -l)
T_NBS=$(echo "$T_NBS / (4 * $PSIZE)" | bc -l)
T_NBE=$(echo "$T_NBE / (4 * $PSIZE)" | bc -l)
T_D2N=$(echo "$T_D2N / (4 * $PSIZE)" | bc -l)
T_0=$(echo "$T_0 / (4 * $PSIZE)" | bc -l)

# Run ms; note:
# -en t i x = subpop i to size x * PSIZE at time t and growth rate to 0
# -es t i p = split subpop i into subpop i and new subpop npop+1, 
# -ej t i j = move all lineages in subpop i to subpop j at time t 
echo "# Command = ms $NSAMP $NREPL -T -t $THETA -r $RHO $NSITE -I $NPOP $NCEU $NYRI $NNEA $NDEN g -en 0 1 1 -es $T_A 1 $AFRAC -ej $T_A 5 3 -ej $T_Y2C 2 1 -en $T_NBS 3 $BFRAC -en $T_NBE 3 1 -ej $T_D2N 4 3 -ej $T_0 3 1"
ms $NSAMP $NREPL -T \
   -t $THETA \
   -r $RHO $NSITE \
   -I $NPOP $NCEU $NYRI $NNEA $NDEN g \
   -en 0 1 1 \
   -es $T_A 1 $AFRAC \
   -ej $T_A 5 3 \
   -ej $T_Y2C 2 1 \
   -en $T_NBS 3 $BFRAC \
   -en $T_NBE 3 1 \
   -ej $T_D2N 4 3 \
   -ej $T_0 3 1 \
   | tail -n1

