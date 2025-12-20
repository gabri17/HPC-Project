#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>
#include <omp.h>

int getProblemSize()
{
    char *iter_str = getenv("ITER");
    if (iter_str == NULL)
    {
        return 10000;
    }
    else
    {
        return atoi(iter_str);
    }
}

int getThreadsNo()
{
    char *iter_str = getenv("OMP_NUM_THREADS");
    if (iter_str == NULL)
    {
        return 1;
    }
    else
    {
        return atoi(iter_str);
    }
}

double f(double x, double y)
{

    int ITER = getProblemSize();

    double dummy = 0;
    int i;
    for (i = 0; i < ITER; i++)
    {
        dummy += 1.0 / (pow(x, 2) + 1.0) + 1.0 / (pow(y, 2) + 1.0);
    }

    // The actual math result: x^5 + y^5
    // We add (dummy * 1e-15) so the compiler does NOT optimize away the loop above.
    double result = pow(x, 5) + pow(y, 5) + (dummy * 1e-15);
    return result;
}

typedef struct
{
    double x, y;
} Point;

/* P = A + u*(B-A) + v*(C-A) */
Point tri_interp(Point A, Point B, Point C, double u, double v)
{
    Point P;
    double xAB = B.x - A.x, yAB = B.y - A.y;
    double xAC = C.x - A.x, yAC = C.y - A.y;

    P.x = A.x + u * xAB + v * xAC;
    P.y = A.y + u * yAB + v * yAC;
    return P;
}

/* Linear interpolation between two points */
Point lerp(Point A, Point B, double t)
{
    Point P;
    P.x = A.x + t * (B.x - A.x);
    P.y = A.y + t * (B.y - A.y);
    return P;
}

/*
  Generate edge nodes (excluding the original 3 vertices) and interior nodes
  for subdivision level nm.
*/
void generate_Em_Im(Point A, Point B, Point C, int nm,
                    Point **E, int *Ecount,
                    Point **I, int *Icount)
{

    if (nm <= 0)
    {
        *E = NULL;
        *Ecount = 0;
        *I = NULL;
        *Icount = 0;
        return;
    }

    int total_pts = (nm + 1) * (nm + 2) / 2;

    // edge points excluding vertices = 3*nm - 3
    int maxE = (nm >= 1) ? (3 * nm - 3) : 0;

    // interior count = total - edgeIncludingVertices
    int maxI = total_pts - 3 * nm;
    if (maxI < 0)
        maxI = 0;

    Point *edge = NULL;
    Point *inter = NULL;

    if (maxE > 0)
    {
        edge = (Point *)malloc(sizeof(Point) * maxE);
    }

    if (maxI > 0)
    {
        inter = (Point *)malloc(sizeof(Point) * maxI);
    }

    int ec = 0, ic = 0;

    /* --- Compute edge nodes --- */
    if (nm >= 2)
    {
        int t;

        // Use a unique name for edge offsets to avoid scope confusion
        int *edge_offsets = malloc(sizeof(int) * (nm));

        edge_offsets[0] = 0;
        for (t = 1; t <= nm - 1; ++t)
        {
            edge_offsets[t] = 3 * (t - 1);
        }
        ec = edge_offsets[nm - 1];
        ec += 3;

#pragma omp parallel for
        for (t = 1; t <= nm - 1; ++t)
        {
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
    // CRITICAL FIX: Only run this if nm > 2.
    // If nm=1, allocating offset array of size 1 and writing to offset[1] crashes the code.
    // If nm=1 or nm=2, there are 0 interior nodes anyway.
    if (nm > 2)
    {
        int *offset = malloc(nm * sizeof(int));

        // This line was causing the crash when nm=1 because offset has size 1 (index 0 only)
        offset[1] = 0;

        int indx = 0;
        for (indx = 2; indx <= nm - 2; indx++)
        {
            offset[indx] = offset[indx - 1] + (nm - 1 - (indx - 1));
        }

        int ia;
        for (ia = 1; ia <= nm - 2; ++ia)
        {
            int ib;
            for (ib = 1; ib <= nm - 1 - ia; ++ib)
            {

                double u = (double)ia / nm;
                double v = (double)ib / nm;

                Point p = tri_interp(A, B, C, u, v);
                inter[ic++] = p;
            }
        }
        free(offset);
    }

    *E = edge;
    *Ecount = ec;
    *I = inter;
    *Icount = ic;
}

// Compute local length for MPI Scatter
int compute_local_length(int processes, int rank, int total_length)
{
    int local_length = total_length / processes;

    if ((total_length % processes) != 0)
    {
        int left = total_length - (processes * local_length);
        if (rank >= (processes - left))
        {
            local_length += 1;
        }
    }
    return local_length;
}

// Compute counts and displacements for MPI Scatter
void compute_counts_and_displs(int *sendcounts, int *displs, int processes, int total_length)
{
    int j;
    for (j = 0; j < processes; j++)
    {
        sendcounts[j] = 0;
        displs[j] = 0;
    }

    int local_length = total_length / processes;
    if ((total_length % processes) != 0)
    {
        int indexOfProcess = processes - 1;
        int left = total_length - (processes * local_length);
        int t;

        for (t = left; t > 0; t--)
        {
            indexOfProcess = processes - left + t - 1;
            sendcounts[indexOfProcess] = 1;
        }
    }

    displs[0] = 0;
    sendcounts[0] = local_length;
    for (j = 1; j < processes; j++)
    {
        sendcounts[j] += local_length;
        displs[j] = sendcounts[j - 1] + displs[j - 1];
    }
}

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

    int processes, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Custom datatype creation
    MPI_Datatype MPI_POINT;
    int elements_in_struct = 2;
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(Point, x);
    offsets[1] = offsetof(Point, y);
    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_create_struct(elements_in_struct, blocklengths, offsets, types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);

    int MaxLevel = 8;
    int paramIndx = 1;
    Point A = {atof(argv[paramIndx++]), atof(argv[paramIndx++])};
    Point B = {
        atof(argv[paramIndx++]),
        atof(argv[paramIndx++]),
    };
    Point C = {
        atof(argv[paramIndx++]),
        atof(argv[paramIndx++]),
    };

    double **R = NULL;
    double startTime = 0;
    if (rank == 0)
    {
        R = malloc(MaxLevel * sizeof(double *));
        size_t i;
        for (i = 0; i < MaxLevel; i++)
        {
            R[i] = malloc(MaxLevel * sizeof(double));
        }
        startTime = MPI_Wtime();
    }

    double area = 0.5 * fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    double sumC = 0;
    if (rank == 0)
    {
        printf("Area is %f\n", area);
        sumC = f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y);
    }

    int m;
    for (m = 0; m < MaxLevel; m++)
    {

        double acc = 0;
        double n_m = (double)(1 << m);
        int triangles = (int)(pow(n_m, 2));

        double sumE = 0, sumI = 0;
        Point *E = NULL, *I = NULL;

        // VLA used here (valid in C99, be careful with very large process counts)
        int sendcounts[processes], displs[processes];
        int Ec = 0, Ic = 0;

        if (rank == 0)
        {
            generate_Em_Im(A, B, C, (int)n_m, &E, &Ec, &I, &Ic);
            printf("m=%d  nm=%f triangles=%d edge_count (excluding vertices)=%d  interior_count=%d\n", m, n_m, triangles, Ec, Ic);

            compute_counts_and_displs(&sendcounts[0], &displs[0], processes, Ec);
        }

        int maxE = (n_m >= 1) ? (3 * n_m - 3) : 0;
        int local_length = compute_local_length(processes, rank, maxE);
        Point *local_E = malloc(sizeof(Point) * local_length);

        // Scatter Edge Nodes
        MPI_Scatterv(E, sendcounts, displs, MPI_POINT, local_E, local_length, MPI_POINT, 0, MPI_COMM_WORLD);

        double local_sum = 0;
#pragma omp parallel for reduction(+ : local_sum)
        for (int indx = 0; indx < local_length; indx++)
        {
            local_sum += f(local_E[indx].x, local_E[indx].y);
        }
        MPI_Reduce(&local_sum, &sumE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        free(local_E);

        if (rank == 0)
        {
            compute_counts_and_displs(&sendcounts[0], &displs[0], processes, Ic);
        }

        int total_pts = (n_m + 1) * (n_m + 2) / 2;
        int maxI = total_pts - 3 * n_m;
        if (maxI < 0)
            maxI = 0;

        local_length = compute_local_length(processes, rank, maxI);
        Point *local_I = malloc(sizeof(Point) * local_length);

        // Scatter Interior Nodes
        MPI_Scatterv(I, sendcounts, displs, MPI_POINT, local_I, local_length, MPI_POINT, 0, MPI_COMM_WORLD);

        local_sum = 0;
#pragma omp parallel for reduction(+ : local_sum)
        for (int indx = 0; indx < local_length; indx++)
        {
            local_sum += f(local_I[indx].x, local_I[indx].y);
        }
        MPI_Reduce(&local_sum, &sumI, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        free(local_I);

        if (rank == 0)
        {
            sumE *= 3.0;
            sumI *= 6.0;
            acc = (area / (3.0 * triangles)) * (sumC + sumE + sumI);
            R[m][0] = acc;

            printf("R[%d][0] = %f\n", m, R[m][0]);
            free(E);
            free(I);
            printf("\n");
        }
    }

    if (rank == 0)
    {
        // Richardson Extrapolation
        printf("R[0]: %f\n", R[0][0]);
        int m;
        for (m = 1; m < MaxLevel; m++)
        {
            printf("R[%d]: %f ", m, R[m][0]);
            int k;
            for (k = 1; k <= m; k++)
            {
                R[m][k] = R[m][k - 1] + (-R[m - 1][k - 1] + R[m][k - 1]) / (pow(2, k) - 1);
                printf("%f ", R[m][k]);
            }
            printf("\n");
        }

        size_t i;
        for (i = 0; i < MaxLevel; i++)
        {
            free(R[i]);
        }
        free(R);

        double endTime = MPI_Wtime();
        printf("Time elapsed %f, %d processors, %d threads per processor, complexity %d, A(%f, %f), B(%f, %f), C(%f, %f)\n",
               endTime - startTime, processes, getThreadsNo(), getProblemSize(), A.x, A.y, B.x, B.y, C.x, C.y);
    }

    MPI_Type_free(&MPI_POINT);
    MPI_Finalize();

    return 0;
}