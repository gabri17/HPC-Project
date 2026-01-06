"""
Given a file
processes time
... ...

We compute speedup and efficiency.
We plot speedup vs #processes and efficiency vs #processes

In cluster load python-3.10.14 and do pip3 install of what needed
"""
import matplotlib.pyplot as plt
import numpy as np
import sys

# ---------- Check arguments ----------
if len(sys.argv) != 3:
    print(f"Usage: python {sys.argv[0]} <datafile> <size>")
    sys.exit(1)

filename = sys.argv[1]
size = sys.argv[2]

# ---------- Load data ----------
# File format:
# processes time

# Automatically skip header
data = np.loadtxt(filename, skiprows=1)  # Always drops first line

processes = data[:, 0].astype(int)
times = data[:, 1]

# ---------- Compute speedup & efficiency ----------
T1 = times[processes == 1][0]

speedup = T1 / times
efficiency = speedup / processes

# ---------- Ideal curves ----------
ideal_speedup = processes
ideal_efficiency = np.ones_like(processes, dtype=float)

# ---------- Plot Speedup ----------
plt.figure()
plt.plot(processes, speedup, marker='o', label='Measured')
plt.plot(processes, ideal_speedup, linestyle='--', label='Ideal')
plt.xlabel("Number of processes")
plt.ylabel("Speedup")
plt.title("Speedup vs Number of Processes")
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig("speedup.png")

# ---------- Plot Efficiency ----------
plt.figure()
plt.plot(processes, efficiency, marker='o', label='Measured')
plt.plot(processes, ideal_efficiency, linestyle='--', label='Ideal')
plt.axhline(0.7, linestyle=':', label='Acceptable (0.7)')
plt.xlabel("Number of processes")
plt.ylabel("Efficiency")
plt.title("Efficiency vs Number of Processes")
plt.legend()
plt.grid(True)
#plt.show()
plt.savefig("efficiency.png")

output_file = "table.txt"
with open(output_file, "w") as f:
    f.write(f"Size: {size}\n")
    print(f"Size: {size}\n")
    f.write("p\tTime\t\tSpeedup\tEfficiency\n")
    print("p\tTime\t\tSpeedup\tEfficiency")
    for p, t, s, e in zip(processes, times, speedup, efficiency):
        f.write(f"{p}\t{t:.6f}\t{s:.3f}\t{e:.3f}\n")
        print(f"{p}\t{t:.6f}\t{s:.3f}\t{e:.3f}")
