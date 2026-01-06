#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb 
#PBS -l walltime=0:15:00
#PBS -q short_HPC4DS
#PBS -V

# Change to the folder where you submitted the job (where executable is located)
cd $PBS_O_WORKDIR

export ITER=1000000
export OMP_NUM_THREADS=1

module load mpich-3.2

# Run the executable in the current directory
mpiexec -n 1 ./serial_implementation 2 1 4 2.5 6 0.75

# For compilation:
# mpicc -g -Wall -fopenmp -o executable source_code.c -lm -std=c99

# To use other PBS strategies:
# -l place=pack:excl