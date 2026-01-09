#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

// Tag definitions
#define TAG_WORK_U 1
#define TAG_WORK_V 2
#define TAG_STOP 3
#define TAG_RESULT 4

#define MAX_LEVEL 8

// Global complexity variable
int G_COMPLEXITY = 1000000;

typedef struct
{
    double x;
    double y;
} Point;

typedef struct
{
    double sum_edge;
    double sum_interior;
} WorkerResult;

double **allocate_matrix(int n)
{
    double **mat = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
        mat[i] = (double *)calloc(n, sizeof(double));
    return mat;
}

void free_matrix(double **mat, int n)
{
    for (int i = 0; i < n; i++)
        free(mat[i]);
    free(mat);
}

// Heavy function with dynamic complexity
double heavy_f(double x, double y)
{
    double dummy = 0;
    for (int i = 0; i < G_COMPLEXITY; i++)
    {
        dummy += 1.0 / (pow(x, 2) + 1.0) + 1.0 / (pow(y, 2) + 1.0);
    }
    double result = pow(x, 5) + pow(y, 5) + (dummy * 1.0e-15);
    return result;
}

void master_code(int size, int buffer_size, Point v1, Point v2, Point v3)
{
    double **R = allocate_matrix(MAX_LEVEL);

    double *buf_u = (double *)malloc(buffer_size * sizeof(double));
    double *buf_v = (double *)malloc(buffer_size * sizeof(double));
    int count = 0;

    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;

    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;
    int dest_worker = 1;

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);
        count = 0;

        // --- Distribution ---
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                buf_u[count] = (double)i / (double)nm;
                buf_v[count] = (double)j / (double)nm;
                count++;

                if (count == buffer_size)
                {
                    MPI_Send(buf_u, buffer_size, MPI_DOUBLE, dest_worker, TAG_WORK_U, MPI_COMM_WORLD);
                    MPI_Send(buf_v, buffer_size, MPI_DOUBLE, dest_worker, TAG_WORK_V, MPI_COMM_WORLD);

                    count = 0;
                    dest_worker++;
                    if (dest_worker >= size)
                        dest_worker = 1;
                }
            }
        }
        if (count > 0)
        {
            MPI_Send(buf_u, count, MPI_DOUBLE, dest_worker, TAG_WORK_U, MPI_COMM_WORLD);
            MPI_Send(buf_v, count, MPI_DOUBLE, dest_worker, TAG_WORK_V, MPI_COMM_WORLD);
            dest_worker++;
            if (dest_worker >= size)
                dest_worker = 1;
        }

        for (int w = 1; w < size; w++)
            MPI_Send(NULL, 0, MPI_DOUBLE, w, TAG_STOP, MPI_COMM_WORLD);

        for (int w = 1; w < size; w++)
        {
            WorkerResult res;
            MPI_Recv(&res, sizeof(WorkerResult), MPI_BYTE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum_edges += res.sum_edge;
            total_sum_interior += res.sum_interior;
        }

        // Update Table
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        // Extrapolation
        double factor = 4.0;
        for (int k = 1; k <= m; k++)
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
    }

    free(buf_u);
    free(buf_v);
    free_matrix(R, MAX_LEVEL);
}

void worker_code(int buffer_size, Point v1, Point v2, Point v3)
{
    double *buf_u = (double *)malloc(buffer_size * sizeof(double));
    double *buf_v = (double *)malloc(buffer_size * sizeof(double));

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        double local_edge_sum = 0.0;
        double local_int_sum = 0.0;
        MPI_Status status;

        while (1)
        {
            MPI_Recv(buf_u, buffer_size, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == TAG_STOP)
                break;

            int received_count;
            MPI_Get_count(&status, MPI_DOUBLE, &received_count);

            MPI_Recv(buf_v, received_count, MPI_DOUBLE, 0, TAG_WORK_V, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (int k = 0; k < received_count; k++)
            {
                double u = buf_u[k];
                double v = buf_v[k];
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);

                double val = heavy_f(px, py);
                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    local_edge_sum += val;
                else
                    local_int_sum += val;
            }
        }
        WorkerResult res = {local_edge_sum, local_int_sum};
        MPI_Send(&res, sizeof(WorkerResult), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
    free(buf_u);
    free(buf_v);
}

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

int main(int argc, char **argv)
{
    int rank, size;

    int buffer_size = 100;
    if (argc > 1)
        buffer_size = atoi(argv[1]);
    if (argc > 2)
        G_COMPLEXITY = atoi(argv[2]);

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

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

    if (rank == 0)
        printf("BENCH,%d,%d,%d,%f\n", G_COMPLEXITY, buffer_size, size, end - start);

    MPI_Finalize();
    return 0;
}