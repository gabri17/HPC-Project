# Data dependency analysis
Index: automatically managed by the pragma directives.

Dynamic extent = lexical extent + code in subprograms invoked

A reduction variable has both a private  and a shared storage behavior.

WR->FLOW (parallel scan)
RW->ANTI
WW->OUTPUT

# Em generation

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                         |
| u               | 98                        | i             | write  | 98                        | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| pAB-pBC-pCA     | 102-105-108               | i             | write  | 102-105-108               | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| ec              | 103-106-109               | i             | write  | 104-107-110               | i             | read   | no  | FLOW   | NO (within same thread) |
| ec              | 103-106-109               | i             | write  | 103-106-109               | i+1           | read   | yes | FLOW   | YES (depends on prev)   |

edge[ec] no issue because any thread is accessing different memory areas of edge. The issue is when the threads access ec (shared var): because it depends on previous values.

Solution: precompute ec.

# Im generation

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                         |
| u               | 121                       | i             | write  | 121                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| v               | 122                       | i             | write  | 122                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| ib              | 119                       | i             | write  | 119                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| p               | 124                       | i             | write  | 124                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| p               | 124                       | i             | write  | 125                       | i             | read   | no  | FLOW   | NO (within same thread) |
| ic              | 125                       | i             | write  | 125                       | i             | read   | no  | FLOW   | NO (within same thread) |
| ic              | 125                       | i             | write  | 125                       | i+1           | read   | yes | FLOW   | YES (depends on prev)   |

- ic same shit of ec
- inner loop not a problem, executed sequentially on the thread with the variable privates
- ic is a shared variable!

Solution: I am precomputing values for ic. Need another loop, but it's less complex than this one.

# Counts and Displs computation

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                         |
| indexOfProcess  | 176                       | i             | write  | 177                       | i             | read   | no  | FLOW   | NO (within same thread) |
| indexOfProcess  | 176                       | i             | write  | 176                       | i+1           | read   | yes | FLOW   | YES (depends on prev)   |

indexOfProcess as ic. We can elaborate a strategy.
indexOfProcess is processes-1-left+t

# Counts and Displs computation 2

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                           |
| sendcounts[j]   | 187                       | i             | read   | 187                       | i             | write  | no  | ANTI | NO (same thread)            |
| sendcounts[j]   | 187                       | i             | write  | 188                       | i+1           | read   | yes | FLOW | YES                         |
| displs[j-1]     | 188                       | i             | read   | 188                       | i+1           | write  | yes | ANTI | YES                         |

Could split the loop and parallelize the computations of sendcounts[j]. Second remains unparallelized.

# Alloc and dealloc

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                           |

no probs, each thraed accessing different memory areas.

# Serial code of each MPI process
| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                           |
| local_sum       | 282                       | i             | write  | 282                       | i+1           | read   | yes | FLOW   | NO (use REDUCE)           |

f is part of dynamic extent but it's not giving problem since it does not modify local_E[indx] or local_I[indx] (= the points themselves).

Big contribution.

<hr>
<hr>
<hr>

# On master process 

Parallelize for 249? Each thread is computing one? All variables are private so it could be done.


# Richardson Extrapolation
| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |         |                          |
| R[m][.]         | 334                       | i             | write  | 282                       | i+1           | read   | yes | FLOW    | YES              |       |

R[m][k] uses R[m][0], R[m][1], ...., R[m][m-1] (no problem second loop is sequential) AND
R[m-1][0], R[m-1][1], ...., R[m-1][m-1]

k private

This cannot be easily parallelized since we have a flow loop-carried data dependency. We did not spent time in designing solutions in for parallelization of this for (with parallel scan technique) since it's not computational heavy: we just loop on a triangular matrix 8*8 doing summation of previoius cells.