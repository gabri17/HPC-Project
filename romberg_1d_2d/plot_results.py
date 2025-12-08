import matplotlib.pyplot as plt
import numpy as np


processes = [1, 2, 4, 8, 16]

# Time in seconds for each configguration
times_1d = [] 
times_2d = []


speedup_1d = [times_1d[0] / t for t in times_1d]
speedup_2d = [times_2d[0] / t for t in times_2d]

efficiency_1d = [s / p for s, p in zip(speedup_1d, processes)]
efficiency_2d = [s / p for s, p in zip(speedup_2d, processes)]

ideal_speedup = processes  

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.plot(processes, speedup_1d, 'o-', label='1D Romberg', color='blue')
ax1.plot(processes, speedup_2d, 's-', label='2D Triangle', color='green')
ax1.plot(processes, ideal_speedup, '--', label='Ideal Linear', color='gray')
ax1.set_title('Speedup vs Number of Processes')
ax1.set_xlabel('Number of Processes (N)')
ax1.set_ylabel('Speedup')
ax1.set_xticks(processes)
ax1.legend()
ax1.grid(True)

ax2.plot(processes, efficiency_1d, 'o-', label='1D Romberg', color='blue')
ax2.plot(processes, efficiency_2d, 's-', label='2D Triangle', color='green')
ax2.set_title('Efficiency vs Number of Processes')
ax2.set_xlabel('Number of Processes (N)')
ax2.set_ylabel('Efficiency (0-1)')
ax2.set_xticks(processes)
ax2.set_ylim(0, 1.1)
ax2.legend()
ax2.grid(True)

plt.tight_layout()
plt.savefig("romberg_performance.png")
plt.show()

print("Graph saved as 'romberg_performance.png'")