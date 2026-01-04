#!/bin/bash
#PBS -l select=32:ncpus=4:mem=2gb -l place=scatter:excl
#PBS -l walltime=2:30:00
#PBS -q short_HPC4DS
#PBS -V

# Change to the folder where job is submitted (where the executable is located)
cd $PBS_O_WORKDIR

export ITER=8000000
export OMP_NUM_THREADS=4
module load mpich-3.2

# Run the executable in the current directory
mpirun.actual -n 32 ../benchmarking 2 1 4 2.5 6 0.75
#perf stat -e cache-misses,cache-references,cycles,instructions mpiexec -n 1 ./benchmarking 2 1 4 2.5 6 0.75

# For compilation:
# mpicc -Wall -fopenmp -o executable source_code.c -lm -std=c99

# To use other PBS strategies:
# -l place=pack:excl
# -l place=scatter:excl