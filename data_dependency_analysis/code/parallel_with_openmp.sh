#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb -l place=pack:excl
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

# Export environment variables to MPI ranks
#PBS -V

export ITER=10
export OMP_NUM_THREADS=1

module load mpich-3.2
mpiexec -np 1 ./HPC-Project/data_dependency_analysis/code/parallel_with_openmp 2 1 4 2.5 6 0.75