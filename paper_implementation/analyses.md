# Data dependency analysis

**Pattern:** Master-Worker with manual chunking.
**Status:** The Master process contains strict flow dependencies that prevent parallelization of the distribution loop. The Worker process contains standard reduction dependencies.

WR -> FLOW (True Dependency)
RW -> ANTI (Anti Dependency)
WW -> OUTPUT (Output Dependency)

# 1. Master: Work Distribution Loop

**Location:** Nested loops `i, j` inside `master_code`.

| Memory Location | Earlier Statement |           |        | Later Statement |           |        | Loop-carried? | Kind of dataflow | Issue                                             |
| --------------- | ----------------- | --------- | ------ | --------------- | --------- | ------ | ------------- | ---------------- | ------------------------------------------------- |
|                 | Line (approx)     | Iteration | Access | Line (approx)   | Iteration | Access |               |                  |                                                   |
| buffer.u[...]   | 151               | i,j       | write  | 151             | i,j+1     | write  | yes           | OUTPUT           | **NO** (Writes to diff index if count increments) |
| buffer.count    | 151               | i,j       | read   | 153             | i,j       | write  | no            | ANTI             | **NO** (Within same thread)                       |
| buffer.count    | 153               | i,j       | write  | 151             | i,j+1     | read   | yes           | FLOW             | **YES** (Must know prev count to write next u)    |
| dest_worker     | 157               | i,j       | read   | 158             | i,j       | write  | no            | ANTI             | **NO** (Within same thread)                       |
| dest_worker     | 158               | i,j       | write  | 157             | i,j+X     | read   | yes           | FLOW             | **YES** (Round-robin logic depends on prev)       |
| MPI_Send        | 156               | i,j       | call   | 156             | i,j+X     | call   | yes           | FLOW             | **YES** (Socket access serialization)             |

**Conclusion:** The Master loop is **inherently serial**. It is not possible to calculate which worker receives the _next_ point without knowing if the _current_ buffer is full.

# 2. Worker: Computation Loop

**Location:** Inner loop `k` inside `worker_code`.

| Memory Location | Earlier Statement |           |        | Later Statement |           |        | Loop-carried? | Kind of dataflow | Issue                       |
| --------------- | ----------------- | --------- | ------ | --------------- | --------- | ------ | ------------- | ---------------- | --------------------------- |
|                 | Line (approx)     | Iteration | Access | Line (approx)   | Iteration | Access |               |                  |                             |
| u               | 220               | k         | read   | 220             | k+1       | read   | no            | INPUT            | **NO** (Read only)          |
| val             | 225               | k         | write  | 225             | k+1       | write  | yes           | OUTPUT           | **NO** (Private var)        |
| local_edge_sum  | 229               | k         | read   | 229             | k         | write  | no            | ANTI             | **NO** (Same thread)        |
| local_edge_sum  | 229               | k         | write  | 229             | k+1       | read   | yes           | FLOW             | **NO** (Standard Reduction) |

**Conclusion:** This loop is safe to parallelize (e.g., SIMD or OMP) because `local_edge_sum` is a reduction variable.

# 3. Serial: Richardson Extrapolation

**Location:** Final loop `k` inside `master_code`.

| Memory Location | Earlier Statement |           |        | Later Statement |           |        | Loop-carried? | Kind of dataflow | Issue                       |
| --------------- | ----------------- | --------- | ------ | --------------- | --------- | ------ | ------------- | ---------------- | --------------------------- |
|                 | Line              | Iteration | Access | Line            | Iteration | Access |               |                  |                             |
| R[m][k]         | 192               | k         | write  | 192             | k+1       | read   | no            | -                | **NO**                      |
| R[m][k-1]       | 192               | k-1       | write  | 192             | k         | read   | yes           | FLOW             | **YES** (Strict Dependency) |

**Conclusion:** `R[m][k]` depends directly on `R[m][k-1]`. This calculation must remain serial, but the cost is negligible compared to the rest of the program.
