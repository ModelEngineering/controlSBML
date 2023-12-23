#!/bin/bash
# Runs processes in parallel
NUM_PROCESS=10
LAST_MODEL=1200
NUM_PER_PROCESS=$( expr $LAST_MODEL / $NUM_PROCESS )
FIRST_MODEL=1
counter=0
rm /tmp/run_process.out
until [ $counter -gt $LAST_MODEL ];
do
LAST=$( expr $counter + $NUM_PER_PROCESS)
python scripts/process_biomodels.py --first $counter --last $LAST >> /tmp/run_process.out &
counter=$( expr $LAST + 1 )
done
