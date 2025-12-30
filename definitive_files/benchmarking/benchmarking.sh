#!/bin/bash
#PBS -l select=1:ncpus=16:mem=2gb -l place=pack:excl
#PBS -l walltime=0:30:00
#PBS -q short_HPC4DS
#PBS -V

# Change to the folder where job is submitted (where the executable is located)
cd $PBS_O_WORKDIR

export ITER=500000
export OMP_NUM_THREADS=16
export OMP_PLACES=threads # threads run on physical cores
#export OMP_PROC_BIND=TRUE
#export OMP_PLACES=cores
module load mpich-3.2

# Run the executable in the current directory
mpiexec -n 1 ./benchmarking 2 1 4 2.5 6 0.75

# For compilation:
# mpicc -g -Wall -fopenmp -o executable source_code.c -lm -std=c99

# To use other PBS strategies:
# -l place=pack:excl
# -l place=scatter:excl