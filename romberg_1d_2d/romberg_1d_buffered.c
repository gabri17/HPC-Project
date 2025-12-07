#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define BUFFER_SIZE 500
#define TAG_WORK 1
#define TAG_STOP 2

double f(double x)
{
    return sin(x);
}

int main(int argc, char **argv)
{
    int rank, size;
    double a = 0.0, b = M_PI;
    int n = 15;
    double **R = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2)
    {
        printf("Error: This Master-Worker model requires at least 2 processes.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    if (rank == 0)
    {
        R = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++)
            R[i] = (double *)malloc(n * sizeof(double));

        double h = b - a;
        R[0][0] = 0.5 * h * (f(a) + f(b));
        printf("Row 0: %f\n", R[0][0]);

        for (int i = 1; i < n; i++)
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

            printf("Row %d Result: %.10f\n", i, R[i][i]);
        }

        for (int i = 0; i < n; i++)
            free(R[i]);
        free(R);
    }
    else
    {
        while (1)
        {

            break;
        }

        for (int i = 1; i < n; i++)
        {
            double row_partial_sum = 0.0;
            double buffer[BUFFER_SIZE];
            MPI_Status status;

            while (1)
            {
                MPI_Recv(buffer, BUFFER_SIZE, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

                if (status.MPI_TAG == TAG_STOP)
                {
                    break;
                }
                int count;
                MPI_Get_count(&status, MPI_DOUBLE, &count);

                for (int k = 0; k < count; k++)
                {
                    row_partial_sum += f(buffer[k]);
                }
            }

            MPI_Send(&row_partial_sum, 1, MPI_DOUBLE, 0, TAG_RESULT, MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();
    return 0;
}