#!/bin/bash
#PBS -l select=1:ncpus=2:mem=2gb -l place=pack:excl
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

# Export environment variables to MPI ranks
#PBS -V

export ITER=1000000
export OMP_NUM_THREADS=2
export OMP_PLACES=threads

module load openmpi-4.0.4
mpiexec --report-bindings -np 1 --map-by node:pe=2 --bind-to core ./HPC-Project/data_dependency_analysis/code/2/parallel_with_openmp 2 1 4 2.5 6 0.75