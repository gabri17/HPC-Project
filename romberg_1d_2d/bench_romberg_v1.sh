#!/bin/bash
#PBS -l select=1:ncpus=64:mem=4gb
#PBS -l walltime=4:00:00
#PBS -q short_cpuQ
#PBS -N RombergBenchmark
#PBS -V

module load mpich-3.2
cd $PBS_O_WORKDIR


OUTPUT_FILE="benchmark_results_4_5run.csv"
echo "Type,Complexity,BufferSize,Procs,Time" > $OUTPUT_FILE

# --- 1. Strong Scalability ---
# Fixed total problem size
STRONG_COMPLEXITY=1000000

echo "Starting Strong Scalability Tests..."

# Loop Buffer from 5 to 250, stepping by 15
for B in $(seq 5 15 250); do
    for P in 1 2 4 8 16 32 64; do
        echo "Running Strong: Procs=$P, Buf=$B"
        mpirun.actual -n $P ./romberg_bench $B $STRONG_COMPLEXITY | grep "BENCH" | awk -F, -v type="Strong" '{printf "%s,%s,%s,%s,%s\n", type, $2, $3, $4, $5}' >> $OUTPUT_FILE
    done
done

# --- 2. Weak Scalability ---
# Problem size grows with processors
BASE_COMPLEXITY=250000

echo "Starting Weak Scalability Tests..."

# Loop Buffer from 5 to 250, stepping by 15
for B in $(seq 5 15 250); do
    for P in 1 2 4 8 16 32 64; do
        CURRENT_COMPLEXITY=$((BASE_COMPLEXITY * P))
        
        echo "Running Weak: Procs=$P, Buf=$B, Cmplx=$CURRENT_COMPLEXITY"
        mpirun.actual -n $P ./romberg_bench $B $CURRENT_COMPLEXITY | grep "BENCH" | awk -F, -v type="Weak" '{printf "%s,%s,%s,%s,%s\n", type, $2, $3, $4, $5}' >> $OUTPUT_FILE
    done
done

echo "Benchmarking complete. Results saved to $OUTPUT_FILE"
