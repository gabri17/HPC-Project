#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define BUFFER_SIZE 200
#define MAX_LEVEL 15

#define TAG_EDGE_DATA 1
#define TAG_INT_DATA 2
#define TAG_STOP 3
#define TAG_RESULT 4

typedef struct
{
    double x;
    double y;
} Point;
typedef struct
{
    int count;
    double u[BUFFER_SIZE];
    double v[BUFFER_SIZE];
} NodeBuffer;
typedef struct
{
    double sum_edge;
    double sum_interior;
} WorkerResult;

double heavy_f(double x, double y)
{

    double junk_math = 0.0;
    for (int k = 0; k < 100; k++)
    {
        junk_math += sqrt(x * k + 1.0) + log(y * k + 2.0);
    }

    return (sin(x * 100.0) * cos(y * 100.0)) + (junk_math * 1e-15);
}

void master_code(int size, Point v1, Point v2, Point v3)
{
    double **R = (double **)malloc(MAX_LEVEL * sizeof(double *));
    for (int i = 0; i < MAX_LEVEL; i++)
        R[i] = (double *)calloc(MAX_LEVEL, sizeof(double));

    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));

    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;
    printf("Level 0 Result: %.10f\n", R[0][0]);

    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);

        NodeBuffer buffer;
        buffer.count = 0;
        int dest_worker = 1;

        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                buffer.u[buffer.count] = (double)i / (double)nm;
                buffer.v[buffer.count] = (double)j / (double)nm;
                buffer.count++;

                if (buffer.count == BUFFER_SIZE)
                {
                    MPI_Send(&buffer, sizeof(NodeBuffer), MPI_BYTE, dest_worker, TAG_EDGE_DATA, MPI_COMM_WORLD);
                    buffer.count = 0;
                    dest_worker++;
                    if (dest_worker >= size)
                        dest_worker = 1;
                }
            }
        }

        if (buffer.count > 0)
        {
            MPI_Send(&buffer, sizeof(NodeBuffer), MPI_BYTE, dest_worker, TAG_EDGE_DATA, MPI_COMM_WORLD);
            dest_worker++;
            if (dest_worker >= size)
                dest_worker = 1;
        }

        for (int w = 1; w < size; w++)
            MPI_Send(NULL, 0, MPI_BYTE, w, TAG_STOP, MPI_COMM_WORLD);

        for (int w = 1; w < size; w++)
        {
            WorkerResult res;
            MPI_Recv(&res, sizeof(WorkerResult), MPI_BYTE, w, TAG_RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            total_sum_edges += res.sum_edge;
            total_sum_interior += res.sum_interior;
        }

        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        double factor = 4.0;
        for (int k = 1; k <= m; k++)
        {
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
            factor *= 4.0;
        }

        printf("Level %d Result: %.10f\n", m, R[m][m]);
    }

    for (int i = 0; i < MAX_LEVEL; i++)
        free(R[i]);
    free(R);
}

void worker_code(int rank, Point v1, Point v2, Point v3)
{
    for (int m = 1; m < MAX_LEVEL; m++)
    {
        double local_edge_sum = 0.0;
        double local_int_sum = 0.0;
        NodeBuffer buffer;
        MPI_Status status;

        while (1)
        {
            MPI_Recv(&buffer, sizeof(NodeBuffer), MPI_BYTE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            if (status.MPI_TAG == TAG_STOP)
                break;

            for (int k = 0; k < buffer.count; k++)
            {
                double u = buffer.u[k];
                double v = buffer.v[k];
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
        WorkerResult res;
        res.sum_edge = local_edge_sum;
        res.sum_interior = local_int_sum;
        MPI_Send(&res, sizeof(WorkerResult), MPI_BYTE, 0, TAG_RESULT, MPI_COMM_WORLD);
    }
}

void serial_code(Point v1, Point v2, Point v3)
{
    double **R = (double **)malloc(MAX_LEVEL * sizeof(double *));
    for (int i = 0; i < MAX_LEVEL; i++)
        R[i] = (double *)calloc(MAX_LEVEL, sizeof(double));

    double area = 0.5 * fabs((v2.x - v1.x) * (v3.y - v1.y) - (v3.x - v1.x) * (v2.y - v1.y));
    double sum_verts = heavy_f(v1.x, v1.y) + heavy_f(v2.x, v2.y) + heavy_f(v3.x, v3.y);
    R[0][0] = (area / 3.0) * sum_verts;
    printf("Level 0 Result: %.10f\n", R[0][0]);

    double total_sum_edges = 0.0;
    double total_sum_interior = 0.0;

    for (int m = 1; m < MAX_LEVEL; m++)
    {
        long long nm = (1LL << m);
        for (long long i = 0; i <= nm; i++)
        {
            for (long long j = 0; j <= nm - i; j++)
            {
                if ((i == 0 && j == 0) || (i == nm && j == 0) || (i == 0 && j == nm))
                    continue;
                if (i % 2 == 0 && j % 2 == 0)
                    continue;

                double u = (double)i / (double)nm;
                double v = (double)j / (double)nm;
                double px = v1.x + u * (v2.x - v1.x) + v * (v3.x - v1.x);
                double py = v1.y + u * (v2.y - v1.y) + v * (v3.y - v1.y);
                double val = heavy_f(px, py);

                double eps = 1e-9;
                if (u < eps || v < eps || (u + v) > (1.0 - eps))
                    total_sum_edges += val;
                else
                    total_sum_interior += val;
            }
        }
        double term = sum_verts + 3.0 * total_sum_edges + 6.0 * total_sum_interior;
        R[m][0] = (area / (3.0 * (double)(nm * nm))) * term;

        double factor = 4.0;
        for (int k = 1; k <= m; k++)
        {
            R[m][k] = R[m][k - 1] + (R[m][k - 1] - R[m - 1][k - 1]) / (factor - 1.0);
            factor *= 4.0;
        }
        printf("Level %d Result: %.10f\n", m, R[m][m]);
    }

    for (int i = 0; i < MAX_LEVEL; i++)
        free(R[i]);
    free(R);
}

int main(int argc, char **argv)
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Point v1 = {0.0, 0.0};
    Point v2 = {1.0, 0.0};
    Point v3 = {0.0, 1.0};

    double start = MPI_Wtime();

    if (size == 1)
    {
        if (rank == 0)
        {
            printf("Running SERIAL mode (N=1)...\n");
            serial_code(v1, v2, v3);
        }
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
        printf("Processes: %d | Time: %f seconds\n", size, end - start);

    MPI_Finalize();
    return 0;
}