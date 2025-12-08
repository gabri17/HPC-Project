#!/bin/bash
#PBS -l select=2:ncpus=8:mem=2gb
#PBS -l walltime=0:50:00
#PBS -q short_cpuQ
#PBS -N Romberg2D_Perf
#PBS -V

module load mpich-3.2
cd $PBS_O_WORKDIR

for P in 1 2 4 8 16; do    
    mpirun.actual -n $P ./romberg_2d_heavy
    
    echo "" 
done