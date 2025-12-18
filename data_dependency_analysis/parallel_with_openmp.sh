#!/bin/bash
#PBS -l select=4:ncpus=2:mem=2gb
#PBS -l walltime=0:10:00
#PBS -q short_cpuQ
#PBS -V

# Change to the folder where you submitted the job (where 'parallel' is located)
cd $PBS_O_WORKDIR

export ITER=1000000
export OMP_NUM_THREADS=2
export OMP_PROC_BIND=TRUE
export OMP_PLACES=cores

module load openmpi-4.0.4

# Run the executable named 'parallel' in the current directory
mpiexec --report-bindings -np 4 --map-by node:pe=2 ./parallel 2 1 4 2.5 6 0.75