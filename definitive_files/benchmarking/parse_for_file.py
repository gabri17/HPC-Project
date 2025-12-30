import os
import re
import argparse

# Regex patterns
filename_pattern = re.compile(r"^benchmarking\.sh\.o\d+")
min_time_pattern = re.compile(r"min_time=([0-9]*\.?[0-9]+)")
mpi_processes_pattern = re.compile(r"mpi_processes=(\d+)")
openmp_threads_pattern = re.compile(r"openmp_threads=(\d+)")

results_mpi = []  # list of (mpi_processes, min_time)
results_openmp = []  # list of (openmp_processes, min_time)

for filename in os.listdir("."):
    if filename_pattern.match(filename):
        print("Matched file:", filename)
        min_time = None
        mpi_processes = None
        openmp_threads = None

        with open(filename, "r") as f:
            for line in f:
                
                if min_time is None:
                    mt = min_time_pattern.search(line)
                    if mt:
                        min_time = float(mt.group(1))

                if mpi_processes is None:
                    mp = mpi_processes_pattern.search(line)
                    if mp:
                        mpi_processes = int(mp.group(1))

                if openmp_threads is None:
                    omp = openmp_threads_pattern.search(line)
                    if omp:
                        openmp_threads = int(omp.group(1))

                if min_time is not None and mpi_processes is not None and openmp_threads is not None:
                    break

        if min_time is not None and mpi_processes is not None:
            results_mpi.append((mpi_processes, min_time))

        if min_time is not None and openmp_threads is not None:
            results_openmp.append((openmp_threads, min_time))

results_mpi.sort(key=lambda x: x[0])
results_openmp.sort(key=lambda x: x[0])

print(results_mpi)
print(results_openmp)

output_file = "results.txt"
print_threads = False

parser = argparse.ArgumentParser(description="Esempio di gestione argomenti")
parser.add_argument('-t', '--threads', action='store_true', default=False, help='Abilita la stampa dei thread')

args = parser.parse_args()

print_threads = args.threads

with open(output_file, "w") as f:
    f.write("processes time\n")
    if print_threads:
        for processes, time in results_openmp:
            f.write(f"{processes} {time}\n")
    else:
        for processes, time in results_mpi:
            f.write(f"{processes} {time}\n")
