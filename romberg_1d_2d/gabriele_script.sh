#!/bin/bash
#PBS -l select=16:ncpus=8:mem=2gb
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

module load mpich-3.2
mpirun.actual -n 32 ./HPC-Project/romberg_1d_2d/romberg_2d_buffered_heavy