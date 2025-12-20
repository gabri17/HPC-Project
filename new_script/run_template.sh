#!/bin/bash
#PBS -l select=__SELECT__:ncpus=4:mem=2gb
#PBS -l walltime=0:15:00
#PBS -q short_HPC4DS
#PBS -V

cd "$PBS_O_WORKDIR"

module load mpich-3.2

export ITER=1000000
export OMP_NUM_THREADS=4

mpiexec -n __NP__ ./parallel 2 1 4 2.5 6 0.75
