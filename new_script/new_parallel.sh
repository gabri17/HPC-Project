#!/bin/bash
#PBS -l select=32:ncpus=4:mem=2gb 
#PBS -l walltime=0:10:00
#PBS -q short_cpuQ
#PBS -V

export ITER=8000000
export OMP_NUM_THREADS=4

module load mpich-3.2

# Run the executable named 'parallel' in the current directory
mpiexec -n 32 ./HPC-Project/new_script/new_parallel 2 1 4 2.5 6 0.75

# mpicc -g -Wall -fopenmp -o parallel parallel.c -lm -std=c99
# -l place=pack:excl