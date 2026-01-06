# HPC-Project

GitHub repository for HPC Project @ Unitn by @me and @Sano.

Academic year 2025-2026

# Code running

To run the code, relies on the following bash files:
- group_implementation/serial_impl/serial_implementation.sh: to run the SERIAL version
- group_implementation/parallel_impl/parallel_implementation.sh: to run the HYBRID PARALLEL version (MPI + OpenMP)
- group_implementation/benchmarking/benchmarking.sh: to run the SERIAL version

In case of cluster execution, ask for proper number of cores (1 core per every thread).

In any case, export to all MPI processes two environment variables:
- ITER: set the complexity of the integrand function (represented by the number of 'dummy' iterations performed)
- OMP_NUM_THREADS: set the number of threads to generate per each MPI process

Any executable takes 6 double arguments: Ax Ay Bx By Cx Cy, specifically the x and y coordinates for the vertices of the cluster.

To run any of the code with a different integrand function, change the function *double f(double x, double y)* in the source code and compile again, following the instructions in the corresponding bash file.