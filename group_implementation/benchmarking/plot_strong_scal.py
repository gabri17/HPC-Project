import numpy as np
import matplotlib.pyplot as plt

processes = np.array([1, 2, 4, 8, 16, 32, 64], dtype=int)
sizes = np.array([250_000, 500_000, 1_000_000, 2_000_000, 4_000_000, 8_000_000], dtype=int)

# Speedup (S) per size (one row per size, columns = processes)
speedup = np.array([
    [1.00, 1.58, 3.92, 7.66, 12.45, 20.14, 20.07],  # 250000
    [1.00, 1.60, 3.94, 6.29, 13.18, 24.15, 38.70],  # 500000
    [1.00, 1.60, 3.91, 6.73, 14.50, 24.26, 24.00],  # 1000000
    [1.00, 1.97, 3.77, 7.75, 14.48, 23.34, 14.23],  # 2000000
    [1.00, 1.97, 3.13, 4.23, 14.56, 21.12, 29.31],  # 4000000
    [1.00, 1.99, 3.78, 7.65, 14.65, 22.48, 35.99],  # 8000000
], dtype=float)

# Efficiency (E) per size (same ordering)
efficiency = np.array([
    [1.00, 0.79, 0.98, 0.96, 0.78, 0.63, 0.31],
    [1.00, 0.80, 0.98, 0.79, 0.82, 0.76, 0.61],
    [1.00, 0.80, 0.98, 0.84, 0.91, 0.76, 0.38],
    [1.00, 0.99, 0.94, 0.97, 0.91, 0.73, 0.22],
    [1.00, 0.98, 0.78, 0.53, 0.91, 0.66, 0.46],
    [1.00, 0.99, 0.95, 0.96, 0.92, 0.70, 0.56],
], dtype=float)

# Plot: Speedup vs processes (one line per size)
plt.figure(figsize=(9,5))
for i, s in enumerate(sizes):
    plt.plot(processes, speedup[i], marker='o', label=f'{s:,}')
# ideal line y = processes
plt.plot(processes, processes, '--', color='gray', label='Ideal')
#plt.xscale('log', base=2)
plt.xticks(processes, [str(p) for p in processes])
plt.xlabel('Processes')
plt.ylabel('Speedup')
plt.title('Speedup vs Processes')
plt.grid(True)
plt.legend(title='Size', loc='upper left', bbox_to_anchor=(1,1))
plt.tight_layout()
plt.savefig('speedup_vs_processes_strong_scal.png', dpi=150)

# Plot: Efficiency vs processes (one line per size)
plt.figure(figsize=(9,5))
for i, s in enumerate(sizes):
    plt.plot(processes, efficiency[i], marker='o', label=f'{s:,}')
plt.axhline(0.7, color='gray', linestyle='--', label='0.7 threshold')
#plt.xscale('log', base=2)
plt.xticks(processes, [str(p) for p in processes])
plt.xlabel('Processes')
plt.ylabel('Efficiency')
plt.title('Efficiency vs Processes')
plt.ylim(0, 1.05)
plt.grid(True)
plt.legend(title='Size', loc='upper left', bbox_to_anchor=(1,1))
plt.tight_layout()
plt.savefig('efficiency_vs_processes_strong_scal.png', dpi=150)

plt.show()
