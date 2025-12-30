#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define BUFFER_SIZE 5
#define MAX_LEVEL 8

#define TAG_WORK_CHUNK 1
#define TAG_STOP 2
#define TAG_RESULT 3

/*
    Struct defining a point in the space
*/
typedef struct
{
    double x;
    double y;
} Point;

/*
    Struct defining buffer sent to the node. It contains:
    - count: #points sent
    - u: x coordinate of the points sent
    - v: y coordinate of the points sent
*/
typedef struct
{
    int count;
    double u[BUFFER_SIZE];
    double v[BUFFER_SIZE];
} NodeBuffer;

/*
    Struct defining local sum of edge nodes and interior nodes given back by the worker.
*/
typedef struct
{
    double sum_edge;
    double sum_interior;
} WorkerResult;

/*
 * Utility function: it allocates a matrix n*n
*/

double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        mat[i] = (double *)calloc(n, sizeof(double));
    return mat;
}

/*
 * Utility function: it free the memory of a matrix n*n
*/
void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
        free(mat[i]);
    free(mat);
}


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
 *   Function to integrate. Must be consiedered as a black box by our algorithm
*/
double heavy_f(double x, double y)
{
    int ITER = getProblemSize();
    double dummy = 0;

    for (int i = 0; i < ITER; i++)
    {
        dummy += 1.0 / (pow(x, 2) + 1.0) + 1.0 / (pow(y, 2) + 1.0);
    }

    double result = pow(x, 5) + pow(y, 5) + (dummy * 1.0e-15);
    return result;
}

/*
 * Serial implementation of the Romberg integration scheme
 * over a triangular domain.
 * Used as a correctness baseline for the parallel version.
*/
void serial_code(Point v1, Point v2, Point v3)
{
    printf("Running in SERIAL mode (N=1)...\n");
    printf("-------------------------------------------------------------\n");

    /* Romberg extrapolation table */
    double **R = allocate_matrix(MAX_LEVEL);

    /*
     * Level 0:
     * Initial approximation using function values
     * at the three triangle vertices.
    */
    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;
    printf("Row  0 | %10.6f\n", R[0][0]);

    /* Accumulators for edge and interior point contributions */
    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;

    /*
     * Loop over refinement levels.
     * Each level doubles the number of subdivisions
     * along each triangle edge.
    */
    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);

        /*
         * Generate and evaluate all new quadrature points
         * introduced at refinement level m.
        */
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                /* Skip triangle vertices (already included) */
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                
                /* Skip points from previous refinement levels */                
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                
                /*
                 * Barycentric coordinates of the point
                 * and mapping to Cartesian coordinates.
                */
                double u = (double)i / (double)nm;
                double v = (double)j / (double)nm;
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);

                /* Evaluate integrand */
                double val = heavy_f(px, py);

                /*
                 * Classify point as edge or interior.
                 * A tolerance is used to handle floating-point error.
                */
                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    total_sum_edges += val;
                else
                    total_sum_interior += val;
            }
        }

        /* Update Column 0 for this level */
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        /* Richardson Extrapolation */
        double factor = 4.0;
        for (int k = 1; k <= m; k++)
        {
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
            factor *= 4.0;
        }

        /* Print Matrix Row */
        printf("Row %2d | ", m);
        for (int k = 0; k <= m; k++)
            printf("%10.6f ", R[m][k]);
        printf("\n");
    }
    free_matrix(R, MAX_LEVEL);
}

/*
 * Master process (rank 0).
 * Coordinates the parallel computation of the Romberg integration table
 * over a triangular domain using MPI workers.
*/
void master_code(int size, Point v1, Point v2, Point v3) {
    printf("Running in PARALLEL mode (N=%d)...\n", size);
    printf("-------------------------------------------------------------\n");

    /* Romberg extrapolation table */
    double **R = allocate_matrix(MAX_LEVEL);

    /* Compute triangle area */
    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    
    /* Initial approximation: evaluate integrand at vertices only */
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;
    printf("Row  0 | %10.6f\n", R[0][0]);

    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);
        NodeBuffer buffer;
        buffer.count = 0;

        /* Next worker to receive a work chunk (round-robin scheduling) */
        int dest_worker = 1;

        /*
         * Generate all new quadrature points at refinement level m.
         * Both edge and interior nodes are generated here.
         *
         * - Vertex nodes are skipped (already accounted for).
         * - Nodes from previous refinement levels are skipped
         *   (even-even index pairs).
         *
         * Workers are responsible for classifying points as edge or interior and to compute actual coordinates from barycentric ones.
        */
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                /* Skip triangle vertices */
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                
                /* Skip points already included in previous levels */
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                /* Barycentric coordinates of the point */
                buffer.u[buffer.count] = (double)i / (double)nm;
                buffer.v[buffer.count] = (double)j / (double)nm;
                buffer.count++;

                /* Send buffer when full and advance to the next worker */
                if (buffer.count == BUFFER_SIZE)
                {
                    MPI_Send(&buffer, sizeof(NodeBuffer), MPI_BYTE, dest_worker, TAG_WORK_CHUNK, MPI_COMM_WORLD);
                    buffer.count = 0;
                    dest_worker++;
                    if (dest_worker >= size)
                        dest_worker = 1;
                }
            }
        }

        /* Send any remaining points */
        if (buffer.count > 0)
        {
            MPI_Send(&buffer, sizeof(NodeBuffer), MPI_BYTE, dest_worker, TAG_WORK_CHUNK, MPI_COMM_WORLD);
            dest_worker++;
            if (dest_worker >= size)
                dest_worker = 1;
        }

        /*
         * Notify all workers that no more work chunks will be sent
         * for this refinement level.
        */
        for (int w = 1; w < size; w++)
            MPI_Send(NULL, 0, MPI_BYTE, w, TAG_STOP, MPI_COMM_WORLD);

        /*
         * Collect partial sums from workers.
         * Each worker returns the accumulated contribution from
         * edge nodes and interior nodes separately.
        */
        for (int w = 1; w < size; w++)
        {
            WorkerResult res;
            MPI_Recv(&res, sizeof(WorkerResult), MPI_BYTE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum_edges += res.sum_edge;
            total_sum_interior += res.sum_interior;
        }

        /* Update Table */
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        /* Richardson extrapolation */
        for (int k = 1; k <= m; k++)
        {
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (pow(2, k) - 1.0);
        }

        /* Print Matrix Row */
        printf("Row %2d | ", m);
        for (int k = 0; k <= m; k++)
            printf("%10.6f ", R[m][k]);
        printf("\n");
    }
    free_matrix(R, MAX_LEVEL);
}

/*
 * Worker process (rank != 0).
 * Receives batches of quadrature points from the master,
 * evaluates the integrand, and accumulates partial sums
 * for edge and interior contributions at each refinement level.
*/
void worker_code(int rank, Point v1, Point v2, Point v3)
{

    /* Loop over all refinement levels (level 0 handled by master) */
    for (int m = 1; m < MAX_LEVEL; m++)
    {
        /* Local accumulators for this refinement level */
        double local_edge_sum = 0.0;
        double local_int_sum = 0.0;
        NodeBuffer buffer;
        MPI_Status status;

        /*
         * Receive work chunks until a STOP message is received
         * from the master, indicating completion of level m.
        */
        while (1)
        {
            MPI_Recv(&buffer, sizeof(NodeBuffer), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            
            /* End of work for current refinement level */
            if (status.MPI_TAG == TAG_STOP)
                break;

            /*
             * Process all points in the received buffer.
             * Points are given in barycentric coordinates (u, v)
             * relative to triangle (v1, v2, v3).
            */
            for (int k = 0; k < buffer.count; k++) {
                
                double u = buffer.u[k];
                double v = buffer.v[k];

                /* Map barycentric coordinates to Cartesian space */
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);

                /* Evaluate integrand at the mapped point */
                double val = heavy_f(px, py);

                /*
                 * Classify point contribution.
                 * A small tolerance is used to robustly detect
                 * points lying on triangle edges.
                */
                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    local_edge_sum += val;
                else
                    local_int_sum += val;
            }
        }

        /*
         * Send accumulated partial sums back to the master
         * for aggregation and Romberg table update.
        */
        WorkerResult res;
        res.sum_edge = local_edge_sum;
        res.sum_interior = local_int_sum;
        MPI_Send(&res, sizeof(WorkerResult), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv) {
    
    //INITIALIZATION
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int paramIndx = 1;

    if(argc != 7){
        if (rank == 0) {
            fprintf(stderr, "Usage: %s Ax Ay Bx By Cx Cy\n", argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    Point v1 = {
        atof(argv[paramIndx++]), 
        atof(argv[paramIndx++])
    };
    Point v2 = {
        atof(argv[paramIndx++]),
        atof(argv[paramIndx++]),
    };
    Point v3 = {
        atof(argv[paramIndx++]),
        atof(argv[paramIndx++]),
    };

    //All aligned before start
    MPI_Barrier(MPI_COMM_WORLD);

    double start = MPI_Wtime();

    if (size == 1)
    {
        if (rank == 0)
            serial_code(v1, v2, v3);
    }
    else
    {
        if (rank == 0)
            master_code(size, v1, v2, v3);
        else
            worker_code(rank, v1, v2, v3);
    }

    double end = MPI_Wtime();
    double elapsed = end - start;
    if (rank == 0)
        printf("\nProcesses: %d | Time: %f seconds\n", size, elapsed);

    MPI_Finalize();
    return 0;
}