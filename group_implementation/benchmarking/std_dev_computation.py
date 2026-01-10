import os
import re
import numpy as np
import pandas as pd
from collections import defaultdict

directory = "."

# Regex patterns
filename_pattern = re.compile(r"^benchmarking\.sh\.o\d+")
warmup_value = re.compile(r"WARMUP=(\d+)")

data = []

def extract_hpc_results(file_path):
    warmup_numbers = []
    complexity = None
    mpi_procs = None
    in_warmup_zone = False
    
    # Patterns
    complexity_pat = r"complexity\s+([\d.]+)"
    mpi_pat = r"mpi_processes\s*=\s*(\d+)"
    number_pat = r"(?<![a-zA-Z0-9_=])\d+\.\d+(?![a-zA-Z0-9_])" # Targets standalone floats

    with open(file_path, 'r') as f:
        for line in f:
            # 1. Capture Complexity
            if "complexity" in line:
                match = re.search(complexity_pat, line)
                if match: complexity = match.group(1)

            # 2. Capture MPI Processes
            if "mpi_processes" in line:
                match = re.search(mpi_pat, line)
                if match: mpi_procs = match.group(1)

            # 3. Handle Warmup Zone
            if "WARMUP=" in line:
                in_warmup_zone = True
                continue
            if "min_time=" in line:
                in_warmup_zone = False
                continue

            if in_warmup_zone:
                # We extract only the floats in this zone
                matches = re.findall(number_pat, line)
                warmup_numbers.extend(matches)

    return {
        "complexity": complexity,
        "mpi_processes": mpi_procs,
        "timings": [float(x) for x in warmup_numbers]
    }

for root, dirs, files in os.walk(directory):
        for file in files:
            if filename_pattern.match(file):
                directory_name = os.path.basename(root)
                if "MPI" in directory_name:

                    if "scatterexcl" in directory_name:
                        logname = directory_name + '/' + file
                        results = extract_hpc_results(logname)
                        results['scatterexcl'] = True
                        data.append(results)
                    else:
                        logname = directory_name + '/' + file
                        results = extract_hpc_results(logname)
                        results['scatterexcl'] = False
                        data.append(results)

rows = []
for d in data:
    std = np.std(d['timings'])
    if int(d['mpi_processes']) != 64:
        rows.append({
            'Complexity': d['complexity'],
            'MPI_Procs': int(d['mpi_processes']),
            'ScatterExcl': 'Yes' if d['scatterexcl'] else 'No',
            'Std_Dev': std
        })

df = pd.DataFrame(rows)

# Reshape for the requested table style: rows (Comp, MPI) and columns (ScatterExcl Yes/No)
pivot_df = df.pivot_table(index=['Complexity', 'MPI_Procs'], columns='ScatterExcl', values='Std_Dev')
pivot_df = pivot_df.reset_index()

# Save to file
output_path = "std_dev.txt"
with open(output_path, "w") as f:
    f.write("Standard Deviation of Timings\n")
    f.write("="*60 + "\n")
    f.write(pivot_df.to_string(index=False))

print(pivot_df)

# 1. Calculate Standard Deviation for every individual run
results = []
for entry in data:
    std_val = np.std(entry['timings'])
    results.append({
        'complexity': entry['complexity'],
        'scatterexcl': entry['scatterexcl'],
        'std': std_val
    })

# 2. Group by Complexity and ScatterExcl, then average the STDs
df = pd.DataFrame(results)
mean_std_df = df.groupby(['complexity', 'scatterexcl'])['std'].mean().reset_index()

# 3. Pivot for tabular display
pivot_df = mean_std_df.pivot(index='complexity', columns='scatterexcl', values='std')
pivot_df.columns = ['No', 'Yes']
pivot_df = pivot_df.reset_index()

print("TABLE: STANDARD DEVIATION OF UNITED TIMING VECTORS")
# Save to file
output_path = "std_dev_just_complexity.txt"
with open(output_path, "w") as f:
    f.write("Standard Deviation of Timings\n")
    f.write("="*60 + "\n")
    f.write(pivot_df.to_string(index=False))

print(pivot_df)