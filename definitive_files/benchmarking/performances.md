# OpenMP scalability

Increasing threads number from 1 to 16 and PBS script with place=pack:excl.

## Size 250,000

| #processes |    Time    | Speedup | Efficiency |
|------------|------------|---------|------------|
| 1          | 42.348636  | 1.000   | 1.000      |
| 2          | 22.142975  | 1.913   | 0.956      |
| 4          | 11.110235  | 3.812   | 0.953      |
| 6          | 7.535454   | 5.620   | 0.937      |
| 8          | 5.685859   | 7.448   | 0.931      |
| 10         | 4.525930   | 9.357   | 0.936      |
| 12         | 4.063082   | 10.423  | 0.869      |
| 14         | 3.460129   | 12.239  | 0.874      |
| 16         | 3.040866   | 13.927  | 0.870      |

## Size 1,000,000

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 170.991604  | 1.000   | 1.000      |
| 2          | 86.867145   | 1.968   | 0.984      |
| 4          | 44.588662   | 3.835   | 0.959      |
| 6          | 30.199103   | 5.662   | 0.944      |
| 8          | 23.594831   | 7.247   | 0.906      |
| 10         | 19.011949   | 8.994   | 0.899      |
| 12         | 15.837571   | 10.797  | 0.900      |
| 14         | 13.849891   | 12.346  | 0.882      |
| 16         | 12.287443   | 13.916  | 0.870      |

## Observation
Good scalabiltiy capabilities thanks to the fact that we are able to distribute the workload good.
Then not too much degradation also because of cache effects.
Exclusivity to avoid interference and having more stable results (std dev of measurements).

# MPI scalability

We choose 4 threads per MPI process, so we can distribute better the workload of function evaluations across threads,
without having too much overhead and permitting a better MPI scalability (with more processes).

We run both with scatter:excl and without, to see effects of the exclusivity of a nodes (= more stable measurements).
With exclusivity we run only up to 32 MPI processes, because with 64 we should have reserved 64 nodes in an exclusive behavior (= difficult to be scheduled at that point).

## Size 250,000

Without exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 12.510785   | 1.000   | 1.000      |
| 2          | 6.770951    | 1.848   | 0.924      |
| 4          | 3.901911    | 3.206   | 0.802      |
| 8          | 1.863524    | 6.714   | 0.839      |
| 16         | 1.100137    | 11.372  | 0.711      |
| 32         | 0.637970    | 19.610  | 0.613      |
| 64         | 0.623258    | 20.073  | 0.314      |


With exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 10.589289   | 1.000   | 1.000      |
| 2          | 6.719587    | 1.576   | 0.788      |
| 4          | 2.702955    | 3.918   | 0.979      |
| 16         | 0.850649    | 12.448  | 0.778      |
| 32         | 0.525721    | 20.142  | 0.629      |

## Size 500,000

Without exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 23.510670   | 1.000   | 1.000      |
| 2          | 13.940681   | 1.686   | 0.843      |
| 4          | 7.341472    | 3.202   | 0.801      |
| 8          | 3.878416    | 6.062   | 0.758      |
| 16         | 1.885109    | 12.472  | 0.779      |
| 32         | 0.966871    | 24.316  | 0.760      |
| 64         | 0.607478    | 38.702  | 0.605      |

With exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 21.506473   | 1.000   | 1.000      |
| 2          | 13.446666   | 1.599   | 0.800      |
| 4          | 5.464328    | 3.936   | 0.984      |
| 8          | 3.418174    | 6.292   | 0.786      |
| 16         | 1.632277    | 13.176  | 0.823      |
| 32         | 0.890675    | 24.146  | 0.755      |

## Size 1,000,000

Without exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 49.491716   | 1.000   | 1.000      |
| 2          | 27.224500   | 1.818   | 0.909      |
| 4          | 12.828467   | 3.858   | 0.964      |
| 8          | 6.444626    | 7.680   | 0.960      |
| 16         | 3.636664    | 13.609  | 0.851      |
| 32         | 2.681617    | 18.456  | 0.577      |
| 64         | 2.061540    | 24.007  | 0.375      |

With exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 42.949641   | 1.000   | 1.000      |
| 2          | 26.892583   | 1.597   | 0.799      |
| 4          | 10.978125   | 3.912   | 0.978      |
| 8          | 6.382632    | 6.729   | 0.841      |
| 16         | 2.962775    | 14.496  | 0.906      |
| 32         | 1.770660    | 24.256  | 0.758      |

## Size 2,000,000

Without exclusivitiy

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 97.757877   | 1.000   | 1.000      |
| 2          | 57.386489   | 1.703   | 0.852      |
| 4          | 26.683056   | 3.664   | 0.916      |
| 8          | 12.811527   | 7.630   | 0.954      |
| 16         | 7.780652    | 12.564  | 0.785      |
| 32         | 7.247288    | 13.489  | 0.422      |
| 64         | 6.868267    | 14.233  | 0.222      |

With exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 85.975512   | 1.000   | 1.000      |
| 2          | 43.609678   | 1.971   | 0.986      |
| 4          | 22.823493   | 3.767   | 0.942      |
| 8          | 11.092168   | 7.751   | 0.969      |
| 16         | 5.938992    | 14.476  | 0.905      |
| 32         | 3.682951    | 23.344  | 0.730      |

## Size 4,000,000 - run it again with scatter:excl

Without exclusivity
| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 197.580485  | 1.000   | 1.000      |
| 2          | 107.711568  | 1.834   | 0.917      |
| 4          | 57.569215   | 3.432   | 0.858      |
| 8          | 31.675536   | 6.238   | 0.780      |
| 16         | 14.323840   | 13.794  | 0.862      |
| 32         | 10.102576   | 19.557  | 0.611      |
| 64         | 6.742064    | 29.306  | 0.458      |

With exclusivity
| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 169.469610  | 1.000   | 1.000      |
| 2          | 86.191074   | 1.966   | 0.983      |
| 4          | 54.173746   | 3.128   | 0.782      |
| 8          | 40.039508   | 4.233   | 0.529      |
| 16         | 11.639300   | 14.560  | 0.910      |
| 32         | 8.024582    | 21.119  | 0.660      |

# Size 8,000,000

Without exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 444.119123  | 1.000   | 1.000      |
| 2          | 260.577173  | 1.704   | 0.852      |
| 4          | 128.916396  | 3.445   | 0.861      |
| 8          | 61.032061   | 7.277   | 0.910      |
| 16         | 31.695204   | 14.012  | 0.876      |
| 32         | 19.914254   | 22.302  | 0.697      |
| 64         | 12.341586   | 35.986  | 0.562      |

With exclusivity

| #processes |    Time     | Speedup | Efficiency |
|------------|-------------|---------|------------|
| 1          | 340.030995  | 1.000   | 1.000      |
| 2          | 171.147568  | 1.987   | 0.993      |
| 4          | 89.954681   | 3.780   | 0.945      |
| 8          | 44.429740   | 7.653   | 0.957      |
| 16         | 23.215947   | 14.646  | 0.915      |
| 32         | 15.129021   | 22.475  | 0.702      |

## Conclusions

Good strong scalability, because by increasing problem size we improve the performances on higher number of MPI processes.
In any case degradation with around 32 MPI processes.

Interesting to see performances are not so good with 250k and 500k but started to be good with 1M.
Why with 250K and 500K initially bad with 2 processes and then good?