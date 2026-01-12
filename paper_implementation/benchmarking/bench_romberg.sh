#!/bin/bash
#PBS -l select=2:ncpus=64:mem=4gb
#PBS -l walltime=04:00:00
#PBS -q short_cpuQ
#PBS -N RombergOpt2
#PBS -V

module load mpich-3.2
cd $PBS_O_WORKDIR


OUTPUT_FILE="benchmark_opt2_results.csv"
echo "Type,Complexity,BufferSize,Procs,Workers,Time" > $OUTPUT_FILE

# --- CONFIGURATION ---
# We use Workers = 1, 2, 4, 8, 16, 32, 64
# P = Workers + 1. So Max Procs = 65.
WORKER_COUNTS=(1 2 4 8 16 32 64)

# Buffer sizes to sweep
BUFFERS=$(seq 5 15 250)

# Base complexities
STRONG_BASE=1000000
WEAK_BASE_PER_WORKER=250000

# ---------------------------------------------------------
# 1. Strong Scalability (Fixed Total Problem Size)
# ---------------------------------------------------------
echo "Starting Strong Scalability (Option 2)..."

for B in $BUFFERS; do
    for W in "${WORKER_COUNTS[@]}"; do
        # P = Workers + 1 Master
        P=$((W + 1))
        
        echo "Strong: Workers=$W (Procs=$P), Buf=$B"
        
        # We pass 'Workers' ($W) to awk to save it in the CSV
        mpirun.actual -n $P ./romberg_bench $B $STRONG_BASE | grep "BENCH" | awk -F, -v type="Strong" -v w=$W '{printf "%s,%s,%s,%s,%s,%s\n", type, $2, $3, $4, w, $5}' >> $OUTPUT_FILE
    done
done

# ---------------------------------------------------------
# 2. Weak Scalability (Problem Size grows with Workers)
# ---------------------------------------------------------
echo "Starting Weak Scalability (Option 2)..."

for B in $BUFFERS; do
    for W in "${WORKER_COUNTS[@]}"; do
        P=$((W + 1))
        
        # Complexity scales with WORKERS (Actual compute units)
        CURRENT_COMPLEXITY=$((WEAK_BASE_PER_WORKER * W))
        
        echo "Weak: Workers=$W (Procs=$P), Buf=$B, Cmplx=$CURRENT_COMPLEXITY"
        
        mpirun.actual -n $P ./romberg_bench $B $CURRENT_COMPLEXITY | grep "BENCH" | awk -F, -v type="Weak" -v w=$W '{printf "%s,%s,%s,%s,%s,%s\n", type, $2, $3, $4, w, $5}' >> $OUTPUT_FILE
    done
done

echo "Benchmarking complete. Results saved to $OUTPUT_FILE"

# For compilation:
# mpicc -Wall -o executable source_code.c -lm -std=c99
