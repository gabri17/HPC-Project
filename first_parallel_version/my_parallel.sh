#!/bin/bash
#PBS -l select=16:ncpus=8:mem=2gb
# set max execution time
#PBS -l walltime=0:10:00
# set the queue
#PBS -q short_cpuQ
module load mpich-3.2
mpirun.actual -n 16 ./proj/my_parallel