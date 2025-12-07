#!/bin/bash
#PBS -l select=2:ncpus=8:mem=2gb
#PBS -l walltime=0:05:00
#PBS -q short_cpuQ
#PBS -N RombergTest
#PBS -V

module load mpich-3.2
cd $PBS_O_WORKDIR
mpirun.actual -n 16 ./romberg_1D