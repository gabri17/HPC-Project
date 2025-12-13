#!/bin/bash
#PBS -l select=1:ncpus=8:mem=2gb
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

# Export environment variables to MPI ranks
#PBS -V

export ITER=1000000
export OMP_NUM_THREADS=8
export OMP_PLACES=threads

module load mpich-3.2
#mpirun.actual -n 1 ./HPC-Project/data_dependency_analysis/my_parallel_with_openmp 2 1 4 2.5 6 0.75
mpiexec --report-bindings -np 1 --map-by node:pe=8 --bind-to core ./HPC-Project/data_dependency_analysis/my_parallel_with_openmp 2 1 4 2.5 6 0.75
#4 2 8 5 12 1.5