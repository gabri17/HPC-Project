#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// --- Constants & Tags ---
// Tags allow the worker to distinguish what kind of message it's receiving
#define TAG_WORK_U 1 // Message contains the 'u' coordinates
#define TAG_WORK_V 2 // Message contains the 'v' coordinates
#define TAG_STOP 3   // Message telling worker to stop for this level
#define TAG_RESULT 4 // Message containing worker's calculated sums

#define MAX_LEVEL 8 // Depth of the Romberg table (how many times we refine grid)

// Global variable to control how much CPU work 'heavy_f' simulates.
// This allows us to test Strong/Weak scaling via command line args.
int G_COMPLEXITY = 1000000;

// Structure representing a point in 2D space
typedef struct
{
    double x;
    double y;
} Point;

// Structure for Worker to return its local partial sums
typedef struct
{
    double sum_edge;     // Sum of points on the triangle boundary
    double sum_interior; // Sum of points inside the triangle
} WorkerResult;

// --- Helper Functions ---

// Allocates a contiguous 2D array
double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        mat[i] = (double *)calloc(n, sizeof(double));
    return mat;
}

// Frees the 2D array
void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
        free(mat[i]);
    free(mat);
}

// The function to integrate. Includes a loop to simulate high CPU load.
// Complexity is controlled by G_COMPLEXITY.
double heavy_f(double x, double y)
{
    double dummy = 0;
    // Busy-work loop to make this function "heavy"
    for (int i = 0; i < G_COMPLEXITY; i++)
    {
        dummy += 1.0 / (pow(x, 2) + 1.0) + 1.0 / (pow(y, 2) + 1.0);
    }
    double result = pow(x, 5) + pow(y, 5) + (dummy * 1.0e-15);
    return result;
}

// --- MASTER PROCESS ---
// Generates coordinates and distributes them to workers
void master_code(int size, int buffer_size, Point v1, Point v2, Point v3)
{
    double **R = allocate_matrix(MAX_LEVEL);

    // Dynamic allocation of buffers based on command line argument
    double *buf_u = (double *)malloc(buffer_size * sizeof(double));
    double *buf_v = (double *)malloc(buffer_size * sizeof(double));
    int count = 0;

    // Step 0: Calculate initial area and vertices (Level 0)
    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;

    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;
    int dest_worker = 1; // Start Round-Robin with Rank 1

    // Loop through Romberg refinement levels
    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);
        count = 0;

        // --- Distribution Phase ---
        // Generate grid points (u, v) for the current level
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                // Skip vertices and previously calculated points (standard Romberg optimization)
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                // Fill buffers
                buf_u[count] = (double)i / (double)nm;
                buf_v[count] = (double)j / (double)nm;
                count++;

                // If buffers are full, send them to the current worker
                if (count == buffer_size)
                {
                    // Note: We send U and V in separate messages with different tags
                    MPI_Send(buf_u, buffer_size, MPI_DOUBLE, dest_worker, TAG_WORK_U, MPI_COMM_WORLD);
                    MPI_Send(buf_v, buffer_size, MPI_DOUBLE, dest_worker, TAG_WORK_V, MPI_COMM_WORLD);

                    count = 0;
                    dest_worker++;
                    if (dest_worker >= size)
                        dest_worker = 1; // Wrap around to first worker
                }
            }
        }
        // Flush remaining items in the buffer
        if (count > 0)
        {
            MPI_Send(buf_u, count, MPI_DOUBLE, dest_worker, TAG_WORK_U, MPI_COMM_WORLD);
            MPI_Send(buf_v, count, MPI_DOUBLE, dest_worker, TAG_WORK_V, MPI_COMM_WORLD);
            dest_worker++;
            if (dest_worker >= size)
                dest_worker = 1;
        }

        // Send STOP signal to all workers for this level
        for (int w = 1; w < size; w++)
            MPI_Send(NULL, 0, MPI_DOUBLE, w, TAG_STOP, MPI_COMM_WORLD);

        // Collect results from all workers
        for (int w = 1; w < size; w++)
        {
            WorkerResult res;
            MPI_Recv(&res, sizeof(WorkerResult), MPI_BYTE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum_edges += res.sum_edge;
            total_sum_interior += res.sum_interior;
        }

        // --- Math Update Phase ---
        // Update Romberg table with collected sums
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        // Richardson Extrapolation
        double factor = 4.0;
        for (int k = 1; k <= m; k++)
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
    }

    free(buf_u);
    free(buf_v);
    free_matrix(R, MAX_LEVEL);
}

// --- WORKER PROCESS ---
// Receives coordinates, calculates heavy_f, returns sums
void worker_code(int buffer_size, Point v1, Point v2, Point v3)
{
    // Allocate buffers to receive incoming data
    double *buf_u = (double *)malloc(buffer_size * sizeof(double));
    double *buf_v = (double *)malloc(buffer_size * sizeof(double));

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        double local_edge_sum = 0.0;
        double local_int_sum = 0.0;
        MPI_Status status;

        while (1)
        {
            // First, wait for the 'U' coordinates (or a STOP signal)
            MPI_Recv(buf_u, buffer_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == TAG_STOP)
                break;

            // Determine how many items were actually sent
            int received_count;
            MPI_Get_count(&status, MPI_DOUBLE, &received_count);

            // Now expect the 'V' coordinates matching the 'U's
            MPI_Recv(buf_v, received_count, MPI_DOUBLE, 0, TAG_WORK_V, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            // Process the batch
            for (int k = 0; k < received_count; k++)
            {
                double u = buf_u[k];
                double v = buf_v[k];

                // Map logical coords (u,v) to physical triangle (px, py)
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);

                double val = heavy_f(px, py);

                // Classify: Edge vs Interior
                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    local_edge_sum += val;
                else
                    local_int_sum += val;
            }
        }

        // Return local results to master
        WorkerResult res = {local_edge_sum, local_int_sum};
        MPI_Send(&res, sizeof(WorkerResult), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
    free(buf_u);
    free(buf_v);
}

// --- SERIAL VERSION ---
// Standard sequential implementation (used when size=1)
void serial_code(Point v1, Point v2, Point v3)
{
    double **R = allocate_matrix(MAX_LEVEL);
    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);
        double edge = 0, interior = 0;
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                double u = (double)i / nm, v = (double)j / nm;
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);
                double val = heavy_f(px, py);
                if (u < 1e-9 || v < 1e-9 || (u + v) > (1.0 - 1e-9))
                    edge += val;
                else
                    interior += val;
            }
        }
        double term = sum_verts + 3.0 * edge + 6.0 * interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;
        for (int k = 1; k <= m; k++)
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (pow(4, k) - 1.0);
    }
    free_matrix(R, MAX_LEVEL);
}

// --- MAIN ---
int main(int argc, char **argv)
{
    int rank, size;

    // Default parameters
    int buffer_size = 100;

    // Parse arguments: ./prog [buffer_size] [complexity]
    if (argc > 1)
        buffer_size = atoi(argv[1]);
    if (argc > 2)
        G_COMPLEXITY = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Triangle Vertices
    Point v1 = {2.0, 1.0};
    Point v2 = {4.0, 2.5};
    Point v3 = {6.0, 0.75};

    double start = MPI_Wtime();

    if (size == 1)
    {
        serial_code(v1, v2, v3);
    }
    else
    {
        if (rank == 0)
            master_code(size, buffer_size, v1, v2, v3);
        else
            worker_code(buffer_size, v1, v2, v3);
    }

    double end = MPI_Wtime();

    // Output formatted specifically for CSV parsing script
    if (rank == 0)
        printf("BENCH,%d,%d,%d,%f\n", G_COMPLEXITY, buffer_size, size, end - start);

    MPI_Finalize();
    return 0;
}