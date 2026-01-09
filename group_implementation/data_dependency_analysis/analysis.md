# Data dependency analysis
Index: automatically managed by the pragma directives.

Dynamic extent = lexical extent + code in subprograms invoked

A reduction variable has both a private and a shared storage behavior.

WR->FLOW
RW->ANTI
WW->OUTPUT

Analysis done on the file *\group_implementation\serial_impl\serial_implementation.c*

## Em generation

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow  |         Issue           |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|-------------------|-------------------------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |               |                   |                         |
| u               | 211                       | i             | write  | 211                       | i + 1         | write  | yes           | OUTPUT            | NO (put private)        |
| pAB-pBC-pCA     | 215-218-221               | i             | write  | 215-218-221               | i + 1         | write  | yes           | OUTPUT            | NO (put private)        |
| pAB-pBC-pCA     | 215-218-221               | i             | write  | 216-219-222               | i + 1         | read   | yes           | FLOW              | NO (within same thread) |
| ec              | 222                       | i             | write  | 216                       | i+1           | read   | yes           | FLOW              | YES (depends on prev it)|

edge[ec] no issue because any thread is accessing different memory areas of edge. The issue is when the threads access ec (shared var): because it depends on previous values.

Solution: precompute ec.

## Im generation

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                         |
| u               | 234                       | i             | write  | 234                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| v               | 235                       | i             | write  | 235                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| ib              | 232                       | ia            | write  | 232                       | ia + 1        | write  | yes | OUTPUT | NO (put private)        |
| p               | 237                       | i             | write  | 237                       | i + 1         | write  | yes | OUTPUT | NO (put private)        |
| p               | 237                       | i             | write  | 238                       | i             | read   | no  | FLOW   | NO (within same thread) |
| ic              | 238                       | i             | write  | 238                       | i+1           | read   | yes | FLOW   | YES (depends on prev it)|

Issues:
- ic same situation of ec
- inner loop not a problem, executed sequentially on the thread with the variable privates
- ic is a shared variable!

Solution: I am precomputing values for ic. Need another loop, but it's less complex than this one.

## Counts and Displs computation

Source code:

```c
int local_length = total_length / processes;

if ((total_length % processes) != 0) {
    int indexOfProcess = processes - 1;
    int left = total_length - (processes * local_length);
    for (; left > 0; left--) {
        sendcounts[indexOfProcess] = 1;             //line 14
        indexOfProcess = indexOfProcess - 1;        //line 15
    }
}
```


| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow |          Issue          |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------------------------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |               |                  |                         |
| indexOfProcess  | 14                       | i             | read   | 15                       | i             | write  | no            | ANTI             | NO (within same thread) |
| indexOfProcess  | 15                       | i             | 14  | 177                       | i+1           | read   | yes           | FLOW             | YES (depends on prev it)|
| indexOfProcess  | 15                       | i             | write  | 15                       | i+1           | write  | yes           | OUTPUT           | YES (depends on prev it)|

indexOfProcess as ic. We can elaborate a strategy (indexOfProcess is processes-1-left+t).

# Counts and Displs computation 2

Source code:
```c
int local_length = total_length / processes;

displs[0] = 0;
sendcounts[0] = local_length;
for (j = 1; j < processes; j++) {
    sendcounts[j] += local_length;                          //line 24
    displs[j] = sendcounts[j - 1] + displs[j - 1];          //line 25
}
```

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                           |
| sendcounts[i]   | 24                       | i             | read   | 24                       | i             | write  | no  | ANTI | NO (within same thread)            |
| sendcounts[i]   | 24                       | i             | write  | 25                       | i+1           | read   | yes | FLOW | YES                         |
| displs[i-1]     | 25                       | i             | read   | 25                       | i+1           | write  | yes | ANTI | YES                         |

Could split the loop and parallelize the computations of sendcounts[j]. Second remains unparallelized.

## Alloc and dealloc

| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                           |

No problems, each thread is accessing different memory areas.

## Serial code of each MPI process
| Memory Location | Earlier Statement         |               |        | Later Statement           |               |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|---------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration     | Access |     |        |                           |
| local_sum       | 336                       | i             | write  | 336                       | i+1           | read   | yes | FLOW   | NO (use REDUCE clause)           |

f is part of dynamic extent but it's not giving problem since it does not modify local_E[indx] or local_I[indx] (= the points themselves).

Here we parallelize using reduction clause.

## Richardson Extrapolation

| Memory Location | Earlier Statement         |               |        | Later Statement           |                  |        | Loop-carried? | Kind of dataflow | Issue |
|-----------------|---------------------------|---------------|--------|---------------------------|------------------|--------|---------------|------------------|-------|
|                 | Line                      | Iteration     | Access | Line                      | Iteration        | Access |     |         |                  |       |
| R[m][k]         | 360                       | m, k          | write  | 360                       | m, k+1           | read   | yes | FLOW    | YES              |       |
| R[m][k]         | 360                       | m, k          | write  | 360                       | m+1, k+1         | read   | yes | FLOW    | YES              |       |

- R[m][k] uses R[m][0], R[m][1], ...., R[m][m-1] (no problem in the second loop if it is sequential) AND
R[m-1][0], R[m-1][1], ...., R[m-1][m-1]
- k private

This cannot be easily parallelized since we have a flow loop-carried data dependency. We did not spent time in designing solutions in for parallelization of this for (with parallel scan technique) since it's not computational heavy: we just loop on a triangular matrix 8x8 doing summation of previoius cells.