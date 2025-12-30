#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>
#include <omp.h>

/*
 *   Utility function: read environment variable ITER and returns the number of iterations needed.
 *  
 *    Returns: 
 *      #iterations of function f (= size of the problem)
*/
int getProblemSize() {
    char *iter_str = getenv("ITER");
    if (iter_str == NULL) {
        return 10000;
    } else {
        return atoi(iter_str);
    }
}

/*
 *   Utility function: read environment variable OMP_NUM_THREADS and returns the number of threads created.
 *  
 *   Returns: 
 *      #threads created per process
*/
int getThreadsNo() {
    char *iter_str = getenv("OMP_NUM_THREADS");
    if (iter_str == NULL) {
        return 1;
    } else {
        return atoi(iter_str);
    }
}

/*
 *   Function to integrate. Must be consiedered as a black box by our algorithm
*/
double f(double x, double y) {

    int ITER = getProblemSize();

    double dummy = 0;
    int i;
    for (i = 0; i < ITER; i++) {
        dummy += 1.0 / (pow(x, 2) + 1.0) + 1.0 / (pow(y, 2) + 1.0);
    }

    // The actual math result: x^5 + y^5
    // We add (dummy * 1e-15) so the compiler does NOT optimize away the loop above.
    double result = pow(x, 5) + pow(y, 5) + (dummy * 1e-15);
    return result;
}

/*
    Struct defining a point in the graph
*/
typedef struct {
    double x, y;
} Point;

/*
 * Perform barycentric (triangular) interpolation inside a triangle.
 *
 * Parameters:
 *   A, B, C : Vertices of the triangle
 *   u, v    : Barycentric interpolation weights
 *
 * Returns:
 *   Interpolated point P
*/
Point tri_interp(Point A, Point B, Point C, double u, double v) {
    /* P = A + u*(B-A) + v*(C-A) 
     * Given triangle vertices A, B, C and barycentric coordinates (u, v),
     * compute the interpolated point P using:
     *
     *     P = A + u * (B - A) + v * (C - A)
     *
     * - When u >= 0, v >= 0, and u + v <= 1, P lies inside the triangle.
     * - u controls movement toward B
     * - v controls movement toward C
    */

    Point P;
    double xAB = B.x - A.x, yAB = B.y - A.y;
    double xAC = C.x - A.x, yAC = C.y - A.y;

    P.x = A.x + u * xAB + v * xAC;
    P.y = A.y + u * yAB + v * yAC;
    return P;
}

/*
 * Perform linear interpolation (LERP) between two points.
 *
 * 
 * Parameters:
 *   A, B : Endpoints of the line segment
 *   t    : Interpolation parameter
 *
 * Returns:
 *   Interpolated point P
 */
Point lerp(Point A, Point B, double t) {
    /*
     * Computes a point P along the line segment from A to B:
     *
     * P = A + t * (B - A)
     *
     * - t = 0 returns A
     * - t = 1 returns B
     * - 0 < t < 1 returns a point between A and B
     * - t < 0 or t > 1 extrapolates beyond the segment
     *
     */

    Point P;
    P.x = A.x + t * (B.x - A.x);
    P.y = A.y + t * (B.y - A.y);
    return P;
}

/*
 * Generate edge and interior nodes for a subdivided triangular region.
 *
 * Given a triangle defined by vertices A, B, and C, this routine generates:
 *  - Edge nodes: points along the edges of the triangle, excluding the
 *    original vertices A, B, and C.
 *  - Interior nodes: points strictly inside the triangle.
 *
 * The density of generated nodes is controlled by the subdivision level nm.
 * Higher values of nm produce finer subdivisions.
 *
 * Memory for the edge node array (*E) and interior node array (*I) is
 * dynamically allocated within this function. The caller is responsible
 * for freeing this memory.
 *
 * Parameters:
 *   A, B, C : Vertices of the triangular region
 *   nm      : Subdivision level
 *   E       : Output pointer to the array of edge nodes
 *   Ecount  : Output count of edge nodes
 *   I       : Output pointer to the array of interior nodes
 *   Icount  : Output count of interior nodes
*/
void generate_Em_Im(Point A, Point B, Point C, int nm,
                    Point **E, int *Ecount,
                    Point **I, int *Icount) {

    if (nm <= 0) {
        *E = NULL;
        *Ecount = 0;
        *I = NULL;
        *Icount = 0;
        return;
    }

    /*
     * Total number of lattice points in a uniformly subdivided triangle.
     *
     * Each edge of the triangle is divided into nm equal segments, producing
     * nm + 1 lattice points along each edge (including endpoints).
     *
     * The resulting triangular lattice contains:
     *
     *     1 + 2 + 3 + ... + (nm + 1)
     *
     * points, corresponding to successive rows of the lattice.
     * This sum evaluates to:
     *
     *     (nm + 1)(nm + 2) / 2
     *
     * which is the total number of lattice points in the subdivided triangle.
     */
    int total_pts = (nm + 1) * (nm + 2) / 2;

    /* edge points including vertices = 3 * nm 
       edge points excluding vertices = 3*nm - 3 (if nm>=1, 0 else) */
    int maxE = (nm >= 1) ? (3 * nm - 3) : 0;
    //each edge gives nm+1 points, so 3(nm+1) - 3 -> each edge include the vertices twices so we count A+...+B+B+...+C+C+...+A -> so we do minus 3
    //3nm - 3, because we do not consider C in Em 

    /* interior count = total - edgeIncludingVertices = total_pts - 3*nm */ //for exclusion
    int maxI = total_pts - 3 * nm;
    if (maxI < 0)
        maxI = 0;

    Point *edge = NULL;
    Point *inter = NULL;

    if (maxE > 0) {
        edge = (Point *)malloc(sizeof(Point) * maxE);
    }

    if (maxI > 0) {
        inter = (Point *)malloc(sizeof(Point) * maxI);
    }

    int ec = 0, ic = 0;

    /* --- Compute edge nodes (exclude the 3 original vertices) --- */
    if (nm >= 2) { /* there are edge points only when nm>=2 */

        int *edge_offsets = malloc(sizeof(int) * (nm));

        edge_offsets[0] = 0;
        for (int t = 1; t <= nm - 1; ++t) {
            edge_offsets[t] = 3 * (t - 1);
        }
        ec = edge_offsets[nm - 1];
        ec += 3;

        #pragma omp parallel for
        for (int t = 1; t <= nm - 1; ++t) {
            double u = (double)t / nm;
            int local_ec = edge_offsets[t];

            Point pAB = lerp(A, B, u);
            edge[local_ec++] = pAB;

            Point pBC = lerp(B, C, u);
            edge[local_ec++] = pBC;

            Point pCA = lerp(C, A, u);
            edge[local_ec++] = pCA;
        }

        free(edge_offsets);
    }

    /* --- Compute interior nodes --- */
    if (nm > 2) { /* there are internal edge points only when nm > 2 */
        int *interior_offsets = malloc(nm * sizeof(int));
        interior_offsets[0] = 0;
        interior_offsets[1] = 0;

        for (int indx = 2; indx <= nm - 2; indx++) {
            interior_offsets[indx] = interior_offsets[indx - 1] + (nm - 1 - (indx - 1));
        }
        ic = interior_offsets[nm - 2] + 1;

        #pragma omp parallel for
        for (int ia = 1; ia <= nm - 2; ++ia) {

            int ic_local = interior_offsets[ia];
            
            for (int ib = 1; ib <= nm - 1 - ia; ++ib) {

                double u = (double)ia / nm;
                double v = (double)ib / nm;

                Point p = tri_interp(A, B, C, u, v);
                inter[ic_local++] = p;
            }
        }
        free(interior_offsets);
    }

    *E = edge;
    *Ecount = ec;
    *I = inter;
    *Icount = ic;
}

/*
 * Compute the number of elements assigned to a specific MPI process.
 *
 * The total dataset of length
 * `total_length` is divided among `processes` processes. This function
 * calculates how many elements the process with rank `rank` should manage.
 *
 * The distribution is as even as possible:
 *  - Each process receives at least floor(total_length / processes) elements.
 *  - If total_length is not divisible by processes, the "leftover" elements
 *    are assigned one by one to the highest-ranked processes.
 *
 * Parameters:
 *   processes    : Total number of MPI processes
 *   rank         : Rank of the current process (0-based)
 *   total_length : Total number of elements to distribute
 *
 * Returns:
 *   The number of elements assigned to the process with the given rank.
*/
int compute_local_length(int processes, int rank, int total_length) {
    int local_length = total_length / processes;

    if ((total_length % processes) != 0) {
        int left = total_length - (processes * local_length);
        if (rank >= (processes - left)) {
            local_length += 1;
        }
    }
    return local_length;
}

/*
 * Compute send counts and displacements for an MPI_Scatterv operation.
 *
 * In MPI, when distributing an array of length `total_length` among
 * `processes` processes, each process may receive a different number of
 * elements if the total length is not divisible evenly. This function
 * computes the two arrays required by MPI_Scatterv:
 *
 *   - sendcounts: number of elements to send to each process
 *   - displs    : displacement (offset) in the source array for each process
 *
 * The distribution is as even as possible:
 *  - Each process receives at least floor(total_length / processes) elements.
 *  - Any remaining "leftover" elements are assigned one by one to the
 *    highest-ranked processes.
 *
 * After computing sendcounts, displacements are calculated cumulatively
 * to indicate the starting index for each process in the original array.
 *
 * Parameters:
 *   sendcounts   : Output array of length `processes` with counts per process
 *   displs       : Output array of length `processes` with displacements
 *   processes    : Total number of MPI processes
 *   total_length : Total number of elements to distribute
 *
 */
void compute_counts_and_displs(int *sendcounts, int *displs, int processes, int total_length) {

    for (int j = 0; j < processes; j++) {
        sendcounts[j] = 0;
        displs[j] = 0;
    }

    int local_length = total_length / processes;
    if ((total_length % processes) != 0) {
        int left = total_length - (processes * local_length);
        int t;

        #pragma omp parallel for
        for (t = left; t > 0; t--) {
            int indexOfProcess = processes - left + t - 1;
            sendcounts[indexOfProcess] = 1;
        }
    }

    displs[0] = 0;
    sendcounts[0] = local_length;
    for (int j = 1; j < processes; j++){
        sendcounts[j] += local_length;
        displs[j] = sendcounts[j - 1] + displs[j - 1];
    }
}

int main(int argc, char *argv[]) {

    //FIRST PART: initialization common to all processes

    MPI_Init(&argc, &argv);

    int processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //we run it WARMUP times, not counting performances
    //then we run it EXECUTION times (taking the minimum execution time)
    int WARMUP = 1;
    int EXECUTION = 10;

    // Custom datatype creation
    MPI_Datatype MPI_POINT;
    int elements_in_struct = 2;
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(Point, x);
    offsets[1] = offsetof(Point, y);
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_create_struct(elements_in_struct, blocklengths, offsets, types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

    int MaxLevel = 8; //maximum level of Richardson extrapolation
    int paramIndx = 1;

    if(argc != 7){
        if (rank == 0) {
            fprintf(stderr, "Usage: %s Ax Ay Bx By Cx Cy\n", argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    Point A = {
        atof(argv[paramIndx++]), 
        atof(argv[paramIndx++])
    };
    Point B = {
        atof(argv[paramIndx++]),
        atof(argv[paramIndx++]),
    };
    Point C = {
        atof(argv[paramIndx++]),
        atof(argv[paramIndx++]),
    };

    if(rank == 0){
        printf("%d processors, %d threads per processor, complexity %d, A(%f, %f), B(%f, %f), C(%f, %f)\n", processes, getThreadsNo(), getProblemSize(), A.x, A.y, B.x, B.y, C.x, C.y);
    }

    int indx = 0;

    double* times = NULL;
    if(rank == 0){
        times = malloc(sizeof(double) * EXECUTION);
    }

    while(indx < (WARMUP + EXECUTION)) {

        //The execution will start at same time
        MPI_Barrier(MPI_COMM_WORLD);

        double **R = NULL;
        double startTime = 0;

        //SECOND PART: initialization of Romberg table by master process and area computation
        if (rank == 0) {
            R = malloc(MaxLevel * sizeof(double *));
            size_t i;
            for (i = 0; i < MaxLevel; i++) {
                R[i] = malloc(MaxLevel * sizeof(double));
            }
            startTime = MPI_Wtime();
        }

        double area = 0.5 * fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
        double sumC = 0;
        if (rank == 0) {
            //printf("Area is %f\n", area);
            sumC = f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y);
        }

        //THIRD PART: computation of first column of Romberg table
        int m;
        for (m = 0; m < MaxLevel; m++) {

            double acc = 0;
            double n_m = (double)(1 << m);
            int triangles = (int)(pow(n_m, 2));

            double sumE = 0, sumI = 0;
            Point *E = NULL, *I = NULL;

            int sendcounts[processes], displs[processes];
            int Ec = 0, Ic = 0;

            //Generation of edge and interior nodes by master process
            //Computation of number and type of edge nodes that each process must manage
            if (rank == 0) {
                generate_Em_Im(A, B, C, (int)n_m, &E, &Ec, &I, &Ic);
                //printf("m=%d  nm=%f triangles=%d edge_count (excluding vertices)=%d  interior_count=%d\n", m, n_m, triangles, Ec, Ic);

                compute_counts_and_displs(&sendcounts[0], &displs[0], processes, Ec);
            }

            //All processes must know the number of edge nodes in total
            int maxE = (n_m >= 1) ? (3 * n_m - 3) : 0;
            int local_length = compute_local_length(processes, rank, maxE);
            Point *local_E = malloc(sizeof(Point) * local_length);

            //Scatter Edge Nodes: distribution of the edge nodes to all the processes 
            MPI_Scatterv(E, sendcounts, displs, MPI_POINT, local_E, local_length, MPI_POINT, 0, MPI_COMM_WORLD);

            //Local evaluation of function f on the edge nodes assigned to the process
            double local_sum = 0;
            #pragma omp parallel for reduction(+ : local_sum)
            for (int indx = 0; indx < local_length; indx++) {
                local_sum += f(local_E[indx].x, local_E[indx].y);
            }

            //Reduction of all local sums to the master process
            MPI_Reduce(&local_sum, &sumE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
            free(local_E);

            //Computation of number and type of interior nodes that each process must manage
            if (rank == 0) {
                compute_counts_and_displs(&sendcounts[0], &displs[0], processes, Ic);
            }

            //All processes must know the number of interior nodes in total
            int total_pts = (n_m + 1) * (n_m + 2) / 2;
            int maxI = total_pts - 3 * n_m;
            if (maxI < 0)
                maxI = 0;

            local_length = compute_local_length(processes, rank, maxI);
            Point *local_I = malloc(sizeof(Point) * local_length);

            //Scatter Interior Nodes: distribution of the interior nodes to all the processes 
            MPI_Scatterv(I, sendcounts, displs, MPI_POINT, local_I, local_length, MPI_POINT, 0, MPI_COMM_WORLD);

            //Local evaluation of function f on the interior nodes assigned to the process
            local_sum = 0;
            #pragma omp parallel for reduction(+ : local_sum)
            for (int indx = 0; indx < local_length; indx++) {
                local_sum += f(local_I[indx].x, local_I[indx].y);
            }

            //Reduction of all local sums to the master process
            MPI_Reduce(&local_sum, &sumI, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            free(local_I);

            //Master process computes and stores the estimation 
            if (rank == 0) {
                sumE *= 3.0;
                sumI *= 6.0;
                acc = (area / (3.0 * triangles)) * (sumC + sumE + sumI);
                R[m][0] = acc;

                //printf("R[%d][0] = %f\n", m, R[m][0]);
                free(E);
                free(I);
                //printf("\n");
            }
        }

        //FOURTH STEP: master process performs richardson extrapolation and computes all remaining entries
        if (rank == 0) {
            
            //Richardson Extrapolation - Print if it's first iteration
            if(indx == 0){
                printf("R[0]: %f\n", R[0][0]);
            }

            int m;
            for (m = 1; m < MaxLevel; m++) {
                if(indx == 0){
                    printf("R[%d]: %f ", m, R[m][0]);
                }
                int k;
                for (k = 1; k <= m; k++) {
                    R[m][k] = R[m][k - 1] + (-R[m - 1][k - 1] + R[m][k - 1]) / (pow(2, k) - 1);
                    if(indx == 0){
                        printf("%f ", R[m][k]);
                    }
                }
                if(indx == 0){
                    printf("\n");
                }
            }

            size_t i;
            for (i = 0; i < MaxLevel; i++) {
                free(R[i]);
            }
            free(R);

            double endTime = MPI_Wtime();

            if(indx >= WARMUP) {
                //We take into account the execution time of the master process as execution time of the entire parallel algorithm
                //We are sure master process will be the last to finish because it must waits for the contribution from all the other processes
                //And then must perform other transformations on the data received.
                times[indx-WARMUP] = endTime - startTime;
            } else {
                printf("\nWARMUP=%f", endTime - startTime);
            }
        }
        indx += 1;
    }
    
    if(rank == 0){
        double minTime = times[0];
        for(int i=0; i < EXECUTION; i++){
            if(times[i] < minTime){
                minTime = times[i];
            }
            printf("\n%f", times[i]);
        }
        printf("\nmin_time=%f", minTime);
        printf("\nmpi_processes=%d", processes);
        printf("\nopenmp_threads=%d", getThreadsNo());

        free(times);
    }


    MPI_Type_free(&MPI_POINT);
    MPI_Finalize();

    return 0;
}