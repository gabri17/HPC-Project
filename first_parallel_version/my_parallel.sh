#!/bin/bash
#PBS -l select=16:ncpus=8:mem=2gb
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ

# Export environment variables to MPI ranks
#PBS -V

export ITER=1000000

module load mpich-3.2
mpirun.actual -n 64 ./proj/my_parallel 4 2 8 5 12 1.5
# 2 1 4 2.5 6 0.75