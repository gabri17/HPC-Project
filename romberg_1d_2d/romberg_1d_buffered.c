#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define BUFFER_SIZE 5
#define TAG_WORK 1
#define TAG_STOP 2
#define TAG_RESULT 3
#define NUM_ROWS 8

// Allocates a square matrix of size n x n
double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        mat[i] = (double *)malloc(n * sizeof(double));
    }
    return mat;
}

// Frees the matrix memory
void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

// 1D function: f(x) = x^5
// Includes a computational delay to simulate "real work"
double heavy_f(double x)
{
    const int ITER = 10000000;
    double dummy = 0;

    // Busy work loop
    for (int i = 0; i < ITER; i++)
    {
        // Use x in the calculation to ensure dependency
        dummy += 1.0 / (pow(x, 2) + 1.0) + sin(i * 0.01);
    }

    // The actual math result: x^5
    // We add (dummy * 1e-15) so the compiler does NOT optimize away the loop above.
    // It is too small to change the actual answer but forces the CPU to do the work.
    double result = pow(x, 5) + (dummy * 1e-15);
    return result;
}

int main(int argc, char **argv)
{
    int rank, size;
    double a = 0.0, b = M_PI; // Integration limits
    double **R = NULL;        // Romberg Table

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double start_time = MPI_Wtime();

    if (size == 1)
    {
        if (rank == 0)
        {
            printf("Running in SERIAL mode (N=1)...\n");
            printf("-------------------------------------------------------------\n");

            R = allocate_matrix(NUM_ROWS);

            // Step 1: Base case (Row 0) - Trapezoidal Rule with 1 segment
            double h = b - a;
            R[0][0] = 0.5 * h * (heavy_f(a) + heavy_f(b));
            printf("Row  0 | %10.6f\n", R[0][0]);

            // Step 2: Iterate through rows
            for (int i = 1; i < NUM_ROWS; i++)
            {
                long long num_new_points = (1LL << (i - 1));
                double h_i = (b - a) / (double)(1LL << i);
                double sum = 0.0;

                // Serial Computation of Function Sums
                for (long long k = 1; k <= num_new_points; k++)
                {
                    double x = a + (2 * k - 1) * h_i;
                    sum += heavy_f(x);
                }

                // Update Column 0
                R[i][0] = 0.5 * R[i - 1][0] + sum * h_i;

                // Richardson Extrapolation
                double factor = 4.0;
                for (int j = 1; j <= i; j++)
                {
                    R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1.0);
                    factor *= 4.0;
                }

                printf("Row %2d | ", i);
                for (int j = 0; j <= i; j++)
                    printf("%10.6f ", R[i][j]);
                printf("\n");
            }
            free_matrix(R, NUM_ROWS);
        }
    }

    else
    {
        if (rank == 0)
        {
            printf("Running in PARALLEL mode (N=%d)...\n", size);
            printf("-------------------------------------------------------------\n");

            R = allocate_matrix(NUM_ROWS);

            // Base case (Calculated locally by Master)
            double h = b - a;
            R[0][0] = 0.5 * h * (heavy_f(a) + heavy_f(b));
            printf("Row  0 | %10.6f\n", R[0][0]);

            // Loop over rows (levels of refinement)
            for (int i = 1; i < NUM_ROWS; i++)
            {
                long long num_new_points = (1LL << (i - 1));
                double h_i = (b - a) / (double)(1LL << i);

                // --- PHASE 1: Distribute Work (Buffered Send) ---
                double buffer[BUFFER_SIZE];
                int count = 0;
                int dest_worker = 1;

                for (long long k = 1; k <= num_new_points; k++)
                {
                    // Generate Coordinate
                    double x = a + (2 * k - 1) * h_i;
                    buffer[count++] = x;

                    // If buffer is full, send to worker
                    if (count == BUFFER_SIZE)
                    {
                        MPI_Send(buffer, BUFFER_SIZE, MPI_DOUBLE, dest_worker, TAG_WORK, MPI_COMM_WORLD);

                        // Round-robin load balancing
                        dest_worker++;
                        if (dest_worker >= size)
                            dest_worker = 1;
                        count = 0;
                    }
                }

                // Send remaining partial buffer
                if (count > 0)
                {
                    MPI_Send(buffer, count, MPI_DOUBLE, dest_worker, TAG_WORK, MPI_COMM_WORLD);
                    dest_worker++; // Keep shifting to keep load balanced
                    if (dest_worker >= size)
                        dest_worker = 1;
                }

                // --- PHASE 2: Signal "End of Row" to Workers ---
                for (int w = 1; w < size; w++)
                {
                    MPI_Send(NULL, 0, MPI_DOUBLE, w, TAG_STOP, MPI_COMM_WORLD);
                }

                // --- PHASE 3: Collect Results (Reduce) ---
                double total_row_sum = 0.0;
                for (int w = 1; w < size; w++)
                {
                    double worker_sum = 0.0;
                    MPI_Recv(&worker_sum, 1, MPI_DOUBLE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    total_row_sum += worker_sum;
                }

                // --- PHASE 4: Update Romberg Table (Serial Part) ---
                R[i][0] = 0.5 * R[i - 1][0] + total_row_sum * h_i;

                // Richardson Extrapolation
                double factor = 4.0;
                for (int j = 1; j <= i; j++)
                {
                    R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (factor - 1.0);
                    factor *= 4.0;
                }

                // Print Row
                printf("Row %2d | ", i);
                for (int j = 0; j <= i; j++)
                    printf("%10.6f ", R[i][j]);
                printf("\n");
            }
            free_matrix(R, NUM_ROWS);
        }
        else
        {
            // Workers loop through rows to stay in sync with Master
            for (int i = 1; i < NUM_ROWS; i++)
            {
                double row_partial_sum = 0.0;
                double buffer[BUFFER_SIZE];
                MPI_Status status;

                // Process buffers until TAG_STOP received
                while (1)
                {
                    MPI_Recv(buffer, BUFFER_SIZE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                    if (status.MPI_TAG == TAG_STOP)
                        break;

                    int count;
                    MPI_Get_count(&status, MPI_DOUBLE, &count);

                    for (int k = 0; k < count; k++)
                    {
                        row_partial_sum += heavy_f(buffer[k]);
                    }
                }
                // Send partial sum for this row back to Master
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