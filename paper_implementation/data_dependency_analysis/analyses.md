# Data dependency analysis

**Pattern:** Master-Worker with manual chunking.
**Status:** The Master process contains strict flow dependencies that prevent parallelization of the distribution loop. The Worker process contains standard reduction dependencies.

WR -> FLOW (True Dependency)
RW -> ANTI (Anti Dependency)
WW -> OUTPUT (Output Dependency)

Analysis done on the file *\paper_implementation\implemenation\romberg_2d_buffered.c*

# 1. Master: Work Distribution Loop

**Location:** Nested loops `i, j` inside `master_code`.

| Memory Location | Earlier Statement |           |        | Later Statement |           |        | Loop-carried? | Kind of dataflow | Issue                                             |
| --------------- | ----------------- | --------- | ------ | --------------- | --------- | ------ | ------------- | ---------------- | ------------------------------------------------- |
|                 | Line (approx)     | Iteration | Access | Line (approx)   | Iteration | Access |               |                  |                                                   |
| buffer.u[...]   | 186               | i,j       | write  | 186             | i,j+1     | write  | yes           | OUTPUT           | **NO** (Writes to diff index if count increments) |
| buffer.v[...]   | 187               | i,j       | write  | 187             | i,j+1     | write  | yes           | OUTPUT           | **NO** (Writes to diff index if count increments) |
| buffer.count    | 186-187           | i,j       | read   | 188             | i,j       | write  | no            | ANTI             | **NO** (Within same thread)                       |
| buffer.count    | 188               | i,j       | write  | 186             | i,j+1     | read   | yes           | FLOW             | **YES** (Must know prev count to write next u)    |
| dest_worker     | 193               | i,j       | read   | 197/199         | i,j       | write  | no            | ANTI             | **NO** (Within same thread)                       |
| dest_worker     | 197/199           | i,j       | write  | 193             | i,j+X     | read   | yes           | FLOW             | **YES** (Round-robin logic depends on prev)       |

**Conclusion:** The Master loop is **inherently serial**. It is not possible to calculate which worker receives the _next_ point without knowing if the _current_ buffer is full.

# 2. Worker: Computation Loop

**Location:** Inner loop `k` inside `worker_code`.

| Memory Location | Earlier Statement |           |        | Later Statement |           |        | Loop-carried? | Kind of dataflow | Issue                       |
| --------------- | ----------------- | --------- | ------ | --------------- | --------- | ------ | ------------- | ---------------- | --------------------------- |
|                 | Line (approx)     | Iteration | Access | Line (approx)   | Iteration | Access |               |                  |                             |
| u               | 271               | k         | write  | 271             | k+1       | write  | yes           | OUTPUT           | **NO** (Private var)        |
| v               | 272               | k         | write  | 272             | k+1       | write  | yes           | OUTPUT           | **NO** (Private var)        |
| val             | 278               | k         | write  | 278             | k+1       | write  | yes           | OUTPUT           | **NO** (Private var)        |
| local_edge_sum  | 283               | k         | read   | 283             | k         | write  | no            | ANTI             | **NO** (Same thread)        |
| local_edge_sum  | 283               | k         | write  | 283             | k+1       | read   | yes           | FLOW             | **NO** (Standard Reduction) |
| local_int_sum   | 285               | k         | read   | 285             | k         | write  | no            | ANTI             | **NO** (Same thread)        |
| local_int_sum   | 285               | k         | write  | 285             | k+1       | read   | yes           | FLOW             | **NO** (Standard Reduction) |

**Conclusion:** This loop is safe to parallelize (e.g., SIMD or OMP) because `local_edge_sum` and `local_int_sum` are reduction variables.

# 3. Serial: Richardson extrapolation

**Location:** Final loop `k` inside `master_code`.

| Memory Location | Earlier Statement |           |        | Later Statement |           |        | Loop-carried? | Kind of dataflow | Issue                       |
| --------------- | ----------------- | --------- | ------ | --------------- | --------- | ------ | ------------- | ---------------- | --------------------------- |
|                 | Line              | Iteration | Access | Line            | Iteration | Access |               |                  |                             |
| R[m][k]         | 234               | k         | write  | 234             | k+1       | read   | yes           | FLOW             | **YES** (Strict Dependency) |

**Conclusion:** `R[m][k]` depends directly on `R[m][k-1]`. This calculation must remain serial, but the cost is negligible compared to the rest of the program.
