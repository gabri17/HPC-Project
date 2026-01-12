#!/bin/bash
#PBS -l select=2:ncpus=8:mem=2gb
#PBS -l walltime=0:20:00
#PBS -q short_cpuQ
#PBS -N Romberg2D_Matrix
#PBS -V

module load mpich-3.2
cd $PBS_O_WORKDIR

mpirun.actual -n 16 ./romberg_2d 2 1 4 2.5 6 0.75

# For compilation:
# mpicc -Wall -o executable source_code.c -lm -std=c99

# To use other PBS strategies:
# -l place=pack:excl
# -l place=scatter:excl