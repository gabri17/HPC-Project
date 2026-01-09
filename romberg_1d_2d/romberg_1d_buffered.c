#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define BUFFER_SIZE 5
#define TAG_WORK 1
#define TAG_STOP 2
#define TAG_RESULT 3
#define NUM_ROWS 8

double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        mat[i] = (double *)malloc(n * sizeof(double));
    }
    return mat;
}

void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
    {
        free(mat[i]);
    }
    free(mat);
}

double heavy_f(double x)
{
    const int ITER = 10000000;
    double dummy = 0;

    // Busy work loop
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
    double a = 0.0, b = M_PI;
    double **R = NULL;

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

            double h = b - a;
            R[0][0] = 0.5 * h * (heavy_f(a) + heavy_f(b));
            printf("Row  0 | %10.6f\n", R[0][0]);

            for (int i = 1; i < NUM_ROWS; i++)
            {
                long long num_new_points = (1LL << (i - 1));
                double h_i = (b - a) / (double)(1LL << i);
                double sum = 0.0;

                for (long long k = 1; k <= num_new_points; k++)
                {
                    double x = a + (2 * k - 1) * h_i;
                    sum += heavy_f(x);
                }

                R[i][0] = 0.5 * R[i - 1][0] + sum * h_i;

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

            double h = b - a;
            R[0][0] = 0.5 * h * (heavy_f(a) + heavy_f(b));
            printf("Row  0 | %10.6f\n", R[0][0]);

            for (int i = 1; i < NUM_ROWS; i++)
            {
                long long num_new_points = (1LL << (i - 1));
                double h_i = (b - a) / (double)(1LL << i);

                double buffer[BUFFER_SIZE];
                int count = 0;
                int dest_worker = 1;

                for (long long k = 1; k <= num_new_points; k++)
                {
                    double x = a + (2 * k - 1) * h_i;
                    buffer[count++] = x;

                    if (count == BUFFER_SIZE)
                    {
                        MPI_Send(buffer, BUFFER_SIZE, MPI_DOUBLE, dest_worker, TAG_WORK, MPI_COMM_WORLD);

                        dest_worker++;
                        if (dest_worker >= size)
                            dest_worker = 1;
                        count = 0;
                    }
                }

                if (count > 0)
                {
                    MPI_Send(buffer, count, MPI_DOUBLE, dest_worker, TAG_WORK, MPI_COMM_WORLD);
                    dest_worker++;
                    if (dest_worker >= size)
                        dest_worker = 1;
                }

                for (int w = 1; w < size; w++)
                {
                    MPI_Send(NULL, 0, MPI_DOUBLE, w, TAG_STOP, MPI_COMM_WORLD);
                }

                double total_row_sum = 0.0;
                for (int w = 1; w < size; w++)
                {
                    double worker_sum = 0.0;
                    MPI_Recv(&worker_sum, 1, MPI_DOUBLE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    total_row_sum += worker_sum;
                }

                R[i][0] = 0.5 * R[i - 1][0] + total_row_sum * h_i;

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
        else
        {
            for (int i = 1; i < NUM_ROWS; i++)
            {
                double row_partial_sum = 0.0;
                double buffer[BUFFER_SIZE];
                MPI_Status status;

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