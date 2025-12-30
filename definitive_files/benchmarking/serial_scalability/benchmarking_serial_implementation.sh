#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb -l place=pack:excl
#PBS -l walltime=4:00:00
#PBS -q short_HPC4DS
#PBS -V

# Change to the folder where job is submitted (where the executable is located)
cd $PBS_O_WORKDIR

export ITER=4000000
export OMP_NUM_THREADS=1
export OMP_PLACES=threads # threads run on physical cores

module load mpich-3.2

# Run the executable in the current directory
mpiexec -n 1 ./benchmarking_serial_implementation 2 1 4 2.5 6 0.75

# For compilation:
# mpicc -g -Wall -fopenmp -o executable source_code.c -lm -std=c99

# To use other PBS strategies:
# -l place=pack:excl
# -l place=scatter:excl