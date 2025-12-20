#!/bin/bash

# List of MPI processes to test
PROCS_LIST=(1 2 4 8 16 32 64)

for NP in "${PROCS_LIST[@]}"; do
    SELECT=$NP  # Number of nodes = number of MPI ranks (adjust if needed)
    
    # Replace placeholders and create a PBS file
    sed -e "s/__SELECT__/$SELECT/" -e "s/__NP__/$NP/" run_template.sh > run_${NP}.sh
    
    # Submit the job
    qsub run_${NP}.sh
done
