#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

// --- Configuration Constants ---
#define BUFFER_SIZE 5 // How many data points to send in one MPI message
#define TAG_WORK 1    // Message tag for work packets
#define TAG_STOP 2    // Message tag to signal end of level
#define TAG_RESULT 3  // Message tag for results
#define NUM_ROWS 8    // Number of Romberg levels

// Helper: Dynamically allocates a 2D array (n x n)
double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        mat[i] = (double *)malloc(n * sizeof(double));
    }
    return mat;
}

// Helper: Frees the memory of the 2D array
void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

// The function we are integrating.
// It includes a busy-wait loop to simulate a computationally expensive task.
double heavy_f(double x)
{
    const int ITER = 10000000;
    double dummy = 0;

    // Simulate heavy CPU load (busy work)
    for (int i = 0; i < ITER; i++)
    {
        dummy += 1.0 / (pow(x, 2) + 1.0) + sin(i * 0.01);
    }

    double result = pow(x, 5) + (dummy * 1e-15);
    return result;
}

int main(int argc, char **argv)
{
    int rank, size;
    double a = 0.0, b = M_PI; // Integration interval [0, PI]
    double **R = NULL;        // The Romberg Table

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();

    // --- CASE 1: Single Processor (Serial Mode) ---
    // If user ran with -n 1, everything runs normally without message passing.
    if (size == 1)
    {
        if (rank == 0)
        {
            printf("Running in SERIAL mode (N=1)...\n");
            printf("-------------------------------------------------------------\n");

            R = allocate_matrix(NUM_ROWS);

            // Step 1: Compute R[0][0] (Trapezoidal rule with 1 interval)
            double h = b - a;
            R[0][0] = 0.5 * h * (heavy_f(a) + heavy_f(b));
            printf("Row  0 | %10.6f\n", R[0][0]);

            // Step 2: Compute subsequent rows
            for (int i = 1; i < NUM_ROWS; i++)
            {
                // Calculate number of new points needed for this level
                long long num_new_points = (1LL << (i - 1));
                double h_i = (b - a) / (double)(1LL << i); // New step size
                double sum = 0.0;

                // Evaluate function at all new midpoints
                for (long long k = 1; k <= num_new_points; k++)
                {
                    double x = a + (2 * k - 1) * h_i;
                    sum += heavy_f(x);
                }

                // Combine with previous result (Trapezoidal update)
                R[i][0] = 0.5 * R[i - 1][0] + sum * h_i;

                // Richardson Extrapolation: Improve accuracy using previous column
                double factor = 4.0;
                for (int j = 1; j <= i; j++)
                {
                    R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1.0);
                    factor *= 4.0;
                }

                // Print the current row of the table
                printf("Row %2d | ", i);
                for (int j = 0; j <= i; j++)
                    printf("%10.6f ", R[i][j]);
                printf("\n");
            }
            free_matrix(R, NUM_ROWS);
        }
    }

    // --- CASE 2: Multiple Processors (Parallel Master-Worker Mode) ---
    else
    {
        // --- MASTER PROCESS (Rank 0) ---
        // Responsible for generating points, sending them to workers, and aggregating results.
        if (rank == 0)
        {
            printf("Running in PARALLEL mode (N=%d)...\n", size);
            printf("-------------------------------------------------------------\n");

            R = allocate_matrix(NUM_ROWS);

            // Compute R[0][0] locally (it's just 2 points, simpler to do here)
            double h = b - a;
            R[0][0] = 0.5 * h * (heavy_f(a) + heavy_f(b));
            printf("Row  0 | %10.6f\n", R[0][0]);

            // Iterate through Romberg levels
            for (int i = 1; i < NUM_ROWS; i++)
            {
                long long num_new_points = (1LL << (i - 1));
                double h_i = (b - a) / (double)(1LL << i);

                // Prepare a buffer to batch multiple points together (reduces network overhead)
                double buffer[BUFFER_SIZE];
                int count = 0;
                int dest_worker = 1; // Start sending to Worker 1

                // Generate points and dispatch to workers
                for (long long k = 1; k <= num_new_points; k++)
                {
                    double x = a + (2 * k - 1) * h_i;
                    buffer[count++] = x;

                    // If buffer is full, send it!
                    if (count == BUFFER_SIZE)
                    {
                        MPI_Send(buffer, BUFFER_SIZE, MPI_DOUBLE, dest_worker, TAG_WORK, MPI_COMM_WORLD);

                        // Round-robin distribution: Next batch goes to the next worker
                        dest_worker++;
                        if (dest_worker >= size)
                            dest_worker = 1; // Wrap around to worker 1 (Rank 0 is Master)
                        count = 0;
                    }
                }

                // Send any remaining points in the partial buffer
                if (count > 0)
                {
                    MPI_Send(buffer, count, MPI_DOUBLE, dest_worker, TAG_WORK, MPI_COMM_WORLD);
                    dest_worker++; // Ensure we don't double-assign if we loop logic later
                    if (dest_worker >= size)
                        dest_worker = 1;
                }

                // Synchronization: Tell ALL workers to stop listening for this level
                for (int w = 1; w < size; w++)
                {
                    MPI_Send(NULL, 0, MPI_DOUBLE, w, TAG_STOP, MPI_COMM_WORLD);
                }

                // Collect results from ALL workers
                double total_row_sum = 0.0;
                for (int w = 1; w < size; w++)
                {
                    double worker_sum = 0.0;
                    // Receive the partial sum calculated by worker 'w'
                    MPI_Recv(&worker_sum, 1, MPI_DOUBLE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    total_row_sum += worker_sum;
                }

                // Update Romberg Table (Trapezoidal Rule)
                R[i][0] = 0.5 * R[i - 1][0] + total_row_sum * h_i;

                // Richardson Extrapolation (Purely serial math, fast enough for Master to do alone)
                double factor = 4.0;
                for (int j = 1; j <= i; j++)
                {
                    R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1.0);
                    factor *= 4.0;
                }

                // Print results
                printf("Row %2d | ", i);
                for (int j = 0; j <= i; j++)
                    printf("%10.6f ", R[i][j]);
                printf("\n");
            }
            free_matrix(R, NUM_ROWS);
        }
        // --- WORKER PROCESSES (Rank > 0) ---
        // Responsible for receiving x-values, computing f(x), and sending back sums.
        else
        {
            // Workers live in a loop for every Level of the Romberg table
            for (int i = 1; i < NUM_ROWS; i++)
            {
                double row_partial_sum = 0.0;
                double buffer[BUFFER_SIZE];
                MPI_Status status;

                // Inner loop: Keep receiving work chunks until Master says STOP
                while (1)
                {
                    // Wait for a message from Master (Rank 0)
                    MPI_Recv(buffer, BUFFER_SIZE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                    // Check the tag: Is it work or a stop signal?
                    if (status.MPI_TAG == TAG_STOP)
                        break;

                    // How many doubles did we actually receive? (Might be less than BUFFER_SIZE)
                    int count;
                    MPI_Get_count(&status, MPI_DOUBLE, &count);

                    // Do the heavy lifting: Compute f(x) for all received points
                    for (int k = 0; k < count; k++)
                    {
                        row_partial_sum += heavy_f(buffer[k]);
                    }
                }
                // Send the total sum for this level back to Master
                MPI_Send(&row_partial_sum, 1, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);
            }
        }
    }

    double end_time = MPI_Wtime();
    if (rank == 0)
    {
        printf("\nTotal Execution Time: %f seconds\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}