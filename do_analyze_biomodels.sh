#!/bin/bash
# Launches analysis of biomodels on different processes.
#  $1 Total number of processes
# Set the NUMPROC
# Environment must be initialized by running activate at the root folder.
TOTALMODEL=1100
if test "$1" = ""
then
    NUMPROC=1
else
    NUMPROC=$1
fi
# Calculate the number of models processed in a batch
NUMMODELS=$(( $TOTALMODEL / $NUMPROC ))
# Start the batches
for ((i=1;i<=$NUMPROC;i++)); do
    if test $i -gt $TOTALMODEL
    then
        break
    fi
    j=$(( $i - 1 ))
    FIRST=$(( $j * $NUMMODELS ))
    FIRST=$(( $FIRST + 1 ))
    LAST=$(( $i * $NUMMODELS ))
    python scripts/process_biomodels.py --first $FIRST --last $LAST &> /tmp/biomodels_${FIRST}_${LAST}.out &
done
