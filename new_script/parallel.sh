#!/bin/bash
#PBS -l select=64:ncpus=4:mem=2gb 
#PBS -l walltime=0:15:00
#PBS -q short_cpuQ
#PBS -V

# Change to the folder where you submitted the job (where 'parallel' is located)
cd $PBS_O_WORKDIR

export ITER=1000000
export OMP_NUM_THREADS=4

module load mpich-3.2

# Run the executable named 'parallel' in the current directory
mpiexec -n 64 ./parallel 2 1 4 2.5 6 0.75

# mpicc -g -Wall -fopenmp -o parallel parallel.c -lm -std=c99
# -l place=pack:excl