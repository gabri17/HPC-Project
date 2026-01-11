# HPC-Project

GitHub repository for HPC Project @ Unitn by @me and @Sano.

Academic year 2025-2026

# Setup

Hardware setup (Unitn CLUSTER HPC2):
- 422,7 TFLOP of theroetical peak performance CPU
- 7674 cores
- 142 calculation nodes
- 10 Gb/s network (some nodes with Infiniband up to 40 Gb/s and some with Omin Path up to 56 Gb/s)

Software setup:
- MPICH version 3.2.
- Python version 3.10.14.
- GCC version 4.8.5.
- Operating system CentOS Linux 7.9
- Scheduler openPBS 19.1

# Code running - OUR algorithm

To run the code, relies on the following bash files:
- group_implementation/serial_impl/serial_implementation.sh: to run the SERIAL version
- group_implementation/parallel_impl/parallel_implementation.sh: to run the HYBRID PARALLEL version (MPI + OpenMP)
- - Specify proper MPI processes to use and threads to use and ask for proper number of chunks and cores per chunk
- group_implementation/benchmarking/benchmarking.sh: to run the BENCHMARKING

In case of cluster execution, ask for proper number of cores (1 core per every thread).

In any case, export to all MPI processes two environment variables:
- ITER: set the complexity of the integrand function (represented by the number of 'dummy' iterations performed)
- OMP_NUM_THREADS: set the number of threads to generate per each MPI process

Any executable takes 6 double arguments: Ax Ay Bx By Cx Cy, specifically the x and y coordinates for the vertices of the cluster.

To run any of the code with a different integrand function, change the function *double f(double x, double y)* in the source code and compile again, following the instructions in the corresponding bash file.

# Code running - Master-Worker algorithm

To run the code, relies on the following bash files:
- paper_implementation/implementation/run_romberg.sh: to run the SERIAL version and the PARALLEL version (MPI only)
- - Specify proper MPI processes to use (1 for serial version) and ask for equal number of cores
- paper_implementation/benchmarking/benchmarking.sh: to run the BENCHMARKING

In case of cluster execution, ask for proper number of cores (1 core per every MPI process).

To run any of the code with a different integrand function, change the function *double f(double x, double y)* in the source code and compile again, following the instructions in the corresponding bash file.
