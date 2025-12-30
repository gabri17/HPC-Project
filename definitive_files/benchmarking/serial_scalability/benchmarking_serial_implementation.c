#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>
#include <omp.h>

/* initial value */
static int GLOBAL_ITER = 1000000;

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

    int ITER = GLOBAL_ITER;

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
    if (nm >= 2) { /* there are internal edge points only when nm>=2 */

        // Use a unique name for edge offsets to avoid scope confusion
        int *edge_offsets = malloc(sizeof(int) * (nm));

        edge_offsets[0] = 0;
        for (int t = 1; t <= nm - 1; ++t) {
            edge_offsets[t] = 3 * (t - 1);
        }
        ec = edge_offsets[nm - 1];
        ec += 3;

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


int main(int argc, char *argv[]) {

    //FIRST PART: initialization

    MPI_Init(&argc, &argv);

    int processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (processes > 1) {
        if(rank == 0){
            fprintf(stderr, "Error! This code must run with only ONE process!");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //we run it WARMUP times, not counting performances
    //then we run it EXECUTION times (taking the minimum execution time)
    int WARMUP = 1;
    int EXECUTION = 10;

    //get size of the problem and store in a global variable
    GLOBAL_ITER = getProblemSize();

    int MaxLevel = 8;
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

    printf("%d processors, %d threads per processor, complexity %d, A(%f, %f), B(%f, %f), C(%f, %f)\n", processes, getThreadsNo(), GLOBAL_ITER, A.x, A.y, B.x, B.y, C.x, C.y);

    int indx = 0;

    double* times = malloc(sizeof(double) * EXECUTION);

    while(indx < (WARMUP + EXECUTION)) {
    
        double **R = NULL;
        double startTime = 0;

        //SECOND PART: initialization of Romberg table and area computation
        R = malloc(MaxLevel * sizeof(double *));
        size_t i;
        for (i = 0; i < MaxLevel; i++) {
            R[i] = malloc(MaxLevel * sizeof(double));
        }
        startTime = MPI_Wtime();
        

        double area = 0.5 * fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
        double sumC = 0;
        //printf("Area is %f\n", area);
        sumC = f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y);

        //THIRD PART: computation of first column of Romberg table
        int m;
        for (m = 0; m < MaxLevel; m++) {

            double acc = 0;
            double n_m = (double)(1 << m);
            int triangles = (int)(pow(n_m, 2));

            double sumE = 0, sumI = 0;
            Point *E = NULL, *I = NULL;

            int Ec = 0, Ic = 0;

            //Generation of edge and interior nodes
            generate_Em_Im(A, B, C, (int)n_m, &E, &Ec, &I, &Ic);
            //printf("m=%d  nm=%f triangles=%d edge_count (excluding vertices)=%d  interior_count=%d\n", m, n_m, triangles, Ec, Ic);

            //Evaluation of function f on the edge nodes
            for (int indx = 0; indx < Ec; indx++) {
                sumE += f(E[indx].x, E[indx].y);
            }

            //Evaluation of function f on the interior nodes
            for (int indx = 0; indx < Ic; indx++) {
                sumI += f(I[indx].x, I[indx].y);
            }

            //Computation of the the estimation 
            sumE *= 3.0;
            sumI *= 6.0;
            acc = (area / (3.0 * triangles)) * (sumC + sumE + sumI);
            R[m][0] = acc;

            //printf("R[%d][0] = %f\n", m, R[m][0]);
            free(E);
            free(I);
            //printf("\n");
        }

        //FOURTH STEP: richardson extrapolation and computation of all remaining entries
            
        //Richardson Extrapolation
        printf("R[0]: %f\n", R[0][0]);
        
        for (m = 1; m < MaxLevel; m++) {
            printf("R[%d]: %f ", m, R[m][0]);
            int k;
            for (k = 1; k <= m; k++) {
                R[m][k] = R[m][k - 1] + (-R[m - 1][k - 1] + R[m][k - 1]) / (pow(2, k) - 1);
                printf("%f ", R[m][k]);
            }
            printf("\n");
        }

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

        indx += 1;
    }

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
    
    MPI_Finalize();

    return 0;
}