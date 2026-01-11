import numpy as np
import matplotlib.pyplot as plt

# Dati forniti
sizes = np.array([250_000, 500_000, 1_000_000, 2_000_000, 4_000_000, 8_000_000], dtype=float)
processes = np.array([1, 2, 4, 8, 16, 32], dtype=float)
times = np.array([10.59, 13.45, 10.98, 11.09, 11.64, 15.13], dtype=float)

# Valori di speedup ed efficienza presi direttamente dalla tabella (non ricalcati)
speedup = np.array([1.00, 1.60, 3.91, 7.75, 14.56, 22.48], dtype=float)
efficiency = np.array([1.00, 0.80, 0.98, 0.97, 0.91, 0.70], dtype=float)

# Plot: Speedup vs processes
plt.figure(figsize=(8,4))
plt.plot(processes, speedup, 'o-', label='Speedup')
# linea ideale y = processes
plt.plot(processes, processes, '--', color='gray', label='Ideal speedup')
plt.xlabel('Processes')
plt.ylabel('Speedup')
plt.title('Speedup vs Processes')
plt.grid(True)
plt.xticks(processes, labels=[str(int(p)) for p in processes])
plt.legend()
plt.tight_layout()
plt.savefig('speedup_vs_processes.png', dpi=150)

# Plot: Efficiency vs processes
plt.figure(figsize=(8,4))
plt.plot(processes, efficiency, 'o-', label='Efficiency')
# linea tratteggiata a 0.7
plt.axhline(0.7, color='gray', linestyle='--', label='0.7 threshold')
plt.xlabel('Processes')
plt.ylabel('Efficiency')
plt.title('Efficiency vs Processes')
plt.grid(True)
plt.xticks(processes, labels=[str(int(p)) for p in processes])
plt.ylim(0,1.05)
plt.legend()
plt.tight_layout()
plt.savefig('efficiency_vs_processes.png', dpi=150)

# Mostra i grafici
plt.show()

# Unica figura con due subplot affiancati
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Speedup
ax = axes[0]
ax.plot(processes, speedup, 'o-', label='Speedup')
ax.plot(processes, processes, '--', color='gray', label='Ideal speedup')
ax.set_xlabel('Processes')
ax.set_ylabel('Speedup')
ax.set_title('Speedup vs Processes')
ax.grid(True)
ax.set_xticks(processes)
ax.set_xticklabels([str(int(p)) for p in processes])
ax.legend()

# Efficiency
ax = axes[1]
ax.plot(processes, efficiency, 'o-', label='Efficiency')
ax.axhline(0.7, color='gray', linestyle='--', label='0.7 threshold')
ax.set_xlabel('Processes')
ax.set_ylabel('Efficiency')
ax.set_title('Efficiency vs Processes')
ax.grid(True)
ax.set_xticks(processes)
ax.set_xticklabels([str(int(p)) for p in processes])
ax.set_ylim(0, 1.05)
ax.legend()

plt.tight_layout()
plt.savefig('weak_scaling_combined.png', dpi=150)
plt.show()
