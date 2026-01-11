#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// --- Configuration ---
#define BUFFER_SIZE 5 // Size of the coordinate buffer sent to workers
#define MAX_LEVEL 8   // Depth of Romberg table

// --- Message Tags ---
#define TAG_WORK_CHUNK 1 // Tag for sending a chunk of work
#define TAG_STOP 2       // Tag to signal end of level
#define TAG_RESULT 3     // Tag for sending back results

// A 2D point structure
typedef struct
{
    double x;
    double y;
} Point;

// A buffer to bundle multiple tasks (coordinates) together
// Sending one by one would be too slow due to network latency.
typedef struct
{
    int count;             // Number of points in this buffer
    double u[BUFFER_SIZE]; // Normalized coordinate u
    double v[BUFFER_SIZE]; // Normalized coordinate v
} NodeBuffer;

// Structure to return results.
// We calculate edges and interior points separately for Romberg logic.
typedef struct
{
    double sum_edge;
    double sum_interior;
} WorkerResult;

// Helper: Allocate 2D array
double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        mat[i] = (double *)calloc(n, sizeof(double));
    return mat;
}

// Helper: Free 2D array
void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
        free(mat[i]);
    free(mat);
}

// The heavy function to integrate f(x,y).
// Includes a loop to simulate computational intensity.
double heavy_f(double x, double y)
{
    const int ITER = 1000000;
    double dummy = 0;

    for (int i = 0; i < ITER; i++)
    {
        dummy += 1.0 / (pow(x, 2) + 1.0) + 1.0 / (pow(y, 2) + 1.0);
    }

    double result = pow(x, 5) + pow(y, 5) + (dummy * 1.0e-15);
    return result;
}

// --- SERIAL VERSION ---
// Runs on a single processor for comparison or if -n 1 is used.
void serial_code(Point v1, Point v2, Point v3)
{
    printf("Running in SERIAL mode (N=1)...\n");
    printf("-------------------------------------------------------------\n");

    double **R = allocate_matrix(MAX_LEVEL);

    // Step 0: Initial Area Calculation (Level 0)
    // Area of triangle using determinant method
    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;
    printf("Row  0 | %10.6f\n", R[0][0]);

    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;

    // Loop through Romberg Levels
    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m); // 2^m segments

        // Generate Grid Points over the unit triangle (0 <= u,v <= 1, u+v <= 1)
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                // Skip the vertices (already calculated)
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                // Optimization: Skip points that were already calculated in previous levels (even indices)
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                // Map normalized (u,v) to physical (x,y) on the triangle
                double u = (double)i / (double)nm;
                double v = (double)j / (double)nm;
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);

                double val = heavy_f(px, py);

                // Classify: Is it on the Edge or Interior?
                // This distinction matters for the specific 2D Romberg formula used
                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    total_sum_edges += val;
                else
                    total_sum_interior += val;
            }
        }

        // Apply 2D Trapezoidal/Romberg Update Formula
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        // Richardson Extrapolation (Refining the result)
        double factor = 4.0;
        for (int k = 1; k <= m; k++)
        {
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
            factor *= 4.0;
        }

        printf("Row %2d | ", m);
        for (int k = 0; k <= m; k++)
            printf("%10.6f ", R[m][k]);
        printf("\n");
    }
    free_matrix(R, MAX_LEVEL);
}

// --- MASTER PROCESS LOGIC ---
// Generates grid points and delegates computation to workers.
void master_code(int size, Point v1, Point v2, Point v3)
{
    printf("Running in PARALLEL mode (N=%d)...\n", size);
    printf("-------------------------------------------------------------\n");

    double **R = allocate_matrix(MAX_LEVEL);

    // Initial Area Calculation (Locally on Master)
    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;
    printf("Row  0 | %10.6f\n", R[0][0]);

    // Accumulators for results across all levels
    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);
        NodeBuffer buffer;
        buffer.count = 0;
        int dest_worker = 1; // Start Round-Robin with Rank 1

        // Loop to generate grid points
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                // Skip vertices and previously calculated points
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                // Pack normalized coordinates into buffer
                buffer.u[buffer.count] = (double)i / (double)nm;
                buffer.v[buffer.count] = (double)j / (double)nm;
                buffer.count++;

                // If buffer is full, ship it to a worker
                if (buffer.count == BUFFER_SIZE)
                {
                    MPI_Send(&buffer, sizeof(NodeBuffer), MPI_BYTE, dest_worker, TAG_WORK_CHUNK, MPI_COMM_WORLD);

                    // Reset buffer and pick next worker
                    buffer.count = 0;
                    dest_worker++;
                    if (dest_worker >= size)
                        dest_worker = 1;
                }
            }
        }
        // Send any remaining items in the buffer
        if (buffer.count > 0)
        {
            MPI_Send(&buffer, sizeof(NodeBuffer), MPI_BYTE, dest_worker, TAG_WORK_CHUNK, MPI_COMM_WORLD);
            dest_worker++;
            if (dest_worker >= size)
                dest_worker = 1;
        }

        // Synchronization: Tell all workers this level is done
        for (int w = 1; w < size; w++)
            MPI_Send(NULL, 0, MPI_BYTE, w, TAG_STOP, MPI_COMM_WORLD);

        // Collect and Aggregate results from all workers
        for (int w = 1; w < size; w++)
        {
            WorkerResult res;
            MPI_Recv(&res, sizeof(WorkerResult), MPI_BYTE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Add worker's partial sums to the global running total
            total_sum_edges += res.sum_edge;
            total_sum_interior += res.sum_interior;
        }

        // Calculate R[m][0] using aggregated sums
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        // Perform Richardson Extrapolation
        double factor = 4.0;
        for (int k = 1; k <= m; k++)
        {
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
            factor *= 4.0;
        }

        // Output results
        printf("Row %2d | ", m);
        for (int k = 0; k <= m; k++)
            printf("%10.6f ", R[m][k]);
        printf("\n");
    }
    free_matrix(R, MAX_LEVEL);
}

// --- WORKER PROCESS LOGIC ---
// Listens for coordinates, computes f(x,y), sums them up.
void worker_code(int rank, Point v1, Point v2, Point v3)
{
    // Workers stay alive for all levels
    for (int m = 1; m < MAX_LEVEL; m++)
    {
        double local_edge_sum = 0.0;
        double local_int_sum = 0.0;
        NodeBuffer buffer;
        MPI_Status status;

        // Processing loop for the current level
        while (1)
        {
            // Wait for data packet
            MPI_Recv(&buffer, sizeof(NodeBuffer), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            // Check if Master signaled end of level
            if (status.MPI_TAG == TAG_STOP)
                break;

            // Process all points in the received buffer
            for (int k = 0; k < buffer.count; k++)
            {
                double u = buffer.u[k];
                double v = buffer.v[k];

                // Map to physical triangle coordinates
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);

                double val = heavy_f(px, py);

                // Classify point (Edge vs Interior)
                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    local_edge_sum += val;
                else
                    local_int_sum += val;
            }
        }

        // Send back the local totals for this level
        WorkerResult res;
        res.sum_edge = local_edge_sum;
        res.sum_interior = local_int_sum;
        MPI_Send(&res, sizeof(WorkerResult), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
}

int main(int argc, char **argv)
{
    int rank, size;

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Define the triangular region
    Point v1 = {2.0, 1.0};
    Point v2 = {4.0, 2.5};
    Point v3 = {6.0, 0.75};

    double start = MPI_Wtime();

    // Select mode based on number of processes
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
    if (rank == 0)
        printf("\nProcesses: %d | Time: %f seconds\n", size, end - start);

    MPI_Finalize();
    return 0;
}