#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include <stddef.h>

#define ITER 800000

double f(double x, double y) {

    double result = 0;
    int i;
    for(i=0; i < ITER; i++){
        result += 1/pow(x, 2) + 1/pow(y, 2) + 1;
    }

    return (result);
}

typedef struct { double x, y; } Point; //struct to represent coordinates of the point

/* P = A + u*(B-A) + v*(C-A) */
Point tri_interp(Point A, Point B, Point C, double u, double v) {
    
    Point P;
    double xAB = B.x - A.x, yAB = B.y - A.y;
    double xAC = C.x - A.x, yAC = C.y - A.y;

    P.x = A.x + u * xAB + v * xAC;
    P.y = A.y + u * yAB + v * yAC;
    return P;
}


/* Linear interpolation between two points */
static Point lerp(Point A, Point B, double t)
{
    Point P;
    P.x = A.x + t * (B.x - A.x); //Px = t Bx + (1-t) Ax
    P.y = A.y + t * (B.y - A.y);
    return P;
}

/*
  Generate edge nodes (excluding the original 3 vertices) and interior nodes
  for subdivision level nm.
  The caller is responsible for freeing *E and *I (malloced here).

  - A, B, C: the vertixes of the triangular region
  - nm: the subdivision level
  - E: the set of edge nodes
  - Ecount: the counts of edge nodes
  - I: the set of interior nodes
  - Icount: the count of interior nodes
*/
void generate_Em_Im(Point A, Point B, Point C, int nm,
                    Point **E, int *Ecount,
                    Point **I, int *Icount) {

    if (nm <= 0) {
        //error: should never happen
        *E = NULL; *Ecount = 0; *I = NULL; *Icount = 0;
        return;
    }

    /* total lattice points = (nm+1)*(nm+2)/2 */
    int total_pts = (nm + 1) * (nm + 2) / 2;
    //each edge is cut into nm equal segments - crucial intuition
    //1 + 2 + 3 + ... + (nm + 1) = (nm+1) * (nm + 2) / 2

    /* edge points including vertices = 3 * nm 
       edge points excluding vertices = 3*nm - 3 (if nm>=1, 0 else) */
    int maxE = (nm >= 1) ? (3 * nm - 3) : 0;
    //each edge gives nm+1 points, so 3(nm+1) - 3 -> each edge include the vertices twices so we count A+...+B+B+...+C+C+...+A -> so we do minus 3
    //3nm - 3, because we do not consider C in Em 

    /* interior count = total - edgeIncludingVertices = total_pts - 3*nm */ //for exclusion
    int maxI = total_pts - 3 * nm;
    if (maxI < 0) maxI = 0; //error: should never happen

    Point *edge = NULL;
    Point *inter = NULL;

    if (maxE > 0) {
        edge = (Point*) malloc(sizeof(Point) * maxE);
    }
    
    if (maxI > 0) {
        inter = (Point*) malloc(sizeof(Point) * maxI);
    }

    int ec = 0, ic = 0;

    /* Compute edge nodes (exclude the 3 original vertices) */
    if (nm >= 2) { /* there are internal edge points only when nm>=2 */
        int t;
        for (t = 1; t <= nm - 1; ++t) { //compute all points on each edge (0...nm) [I have nm+1 points][0 is A so I exclude it][nm is B so I exclude it]
            double u = (double) t / nm; //contribution of each point
            
            //each point on the edge is a linear contribution of its extremes

            Point pAB = lerp(A, B, u);
            edge[ec++] = pAB;

            Point pBC = lerp(B, C, u);
            edge[ec++] = pBC;

            Point pCA = lerp(C, A, u);
            edge[ec++] = pCA;
        }
    }

    /* Interior nodes: ia,ib,ic >=1 and ia+ib+ic = nm */
    int ia;
    for (ia = 1; ia <= nm - 2; ++ia) {
        /* ib can range from 1 to nm-1-ia, such that ic = nm-ia-ib >=1 */
        int ib;
        for (ib = 1; ib <= nm - 1 - ia; ++ib) {
            
            double u = (double) ia / nm;
            double v = (double) ib / nm;

            Point p = tri_interp(A, B, C, u, v);
            inter[ic++] = p;
        }
    }

    *E = edge;
    *Ecount = ec;
    *I = inter;
    *Icount = ic;
}


//auxiliary functions for nodes distribution
int compute_local_length(int processes, int rank, int total_length) {
	int local_length = total_length / processes;

	if((total_length % processes) != 0){
		int left = total_length - (processes*local_length);
		if(rank >= (processes - left)){
			local_length +=1;
		}
	}

	return local_length;
}

void compute_counts_and_displs(int* sendcounts, int* displs, int processes, int total_length){
    int j;
    for(j = 0; j < processes; j++){
        sendcounts[j] = 0;
        displs[j] = 0;
    }

	int local_length = total_length / processes;
	if((total_length % processes) != 0){
		//alcuni ne avranno 1 in + (sicuramente saranno < di processes)
		int indexOfProcess = processes-1;
		int left = total_length - (processes*local_length);
		//printf("Chi ne ha di piu: ");
		while(indexOfProcess >= 0 && left > 0){
			sendcounts[indexOfProcess] = 1;
			//printf("%d ", indexOfProcess);
			indexOfProcess--;
			left--;
		}
		//printf("\n");
	}

	displs[0] = 0;
	sendcounts[0] = local_length;
	//printf("Primo valore: 0-%d ", sendcounts[0] + displs[0] - 1);
	for(j = 1; j < processes; j++){
		sendcounts[j] += local_length;
		displs[j] = sendcounts[j-1] + displs[j-1];
		//processo 0 i primi local_length, 1 i secondi local_length
		//from rank * local_length to (rank+1) * local_length - 1
		//printf("%d-%d ", displs[j], displs[j] + sendcounts[j] - 1);
	}
	//printf("\n");
}

int main(int argc, char* argv[]){

    MPI_Init(NULL, NULL);

	int processes, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //custom datatype creation
    MPI_Datatype MPI_POINT;

	int elements_in_struct = 2;
    int blocklengths[2] = {1, 1};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(Point, x);
    offsets[1] = offsetof(Point, y);
	//offsetof is a C standard macro that tells you how many bytes from the start of a struct a certain field is located.

    MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};

    MPI_Type_create_struct(elements_in_struct, blocklengths, offsets, types, &MPI_POINT);
    MPI_Type_commit(&MPI_POINT);


    
    //paremeters chosen
    int MaxLevel = 8;

    Point A = {2.0, 1.0};
    Point B = {4.0, 2.5};
    Point C = {6.0, 0.75};


    //INPUT PREPARATION
    double** R = NULL;
    double startTime = 0;
    if(rank == 0){
        R = malloc(MaxLevel * sizeof(double*));
        
        size_t i;
        for(i = 0; i < MaxLevel; i++){
            R[i] = malloc(MaxLevel * sizeof(double));
        }
        startTime = MPI_Wtime();
    }
    
    double area = 0.5 * fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    if(rank == 0){
        printf("Area is %f\n", area);
        //R[0][0] =  (area / 3.0) * (f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y));
    }
    
    int m;
    for(m = 0; m < MaxLevel; m++){

        //compute R[0][0], R[1][0], ..., R[m][0], ..., R[MaxLevel][0]
        double acc = 0;
        double n_m = (double) (1 << m);
        int triangles = (int) (pow(n_m, 2));

        double sumC = 0, sumE = 0, sumI = 0; //C, Em, Im
        Point *E = NULL, *I = NULL;
        double *EflattenX = NULL, *EflattenY = NULL;
        double *IflattenX = NULL, *IflattenY = NULL;

        int sendcounts[processes], displs[processes];
        int j;
        int Ec = 0, Ic = 0;
        for(j = 0; j < processes; j++){
            sendcounts[j] = 0;
            displs[j] = 0;
        }


        if(rank == 0){
            //first sum member
            sumC += f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y);

            //edge and interior nodes computation
            generate_Em_Im(A, B, C, (int) n_m, &E, &Ec, &I, &Ic);
            printf("m=%d  nm=%f triangles=%d edge_count (excluding vertices)=%d  interior_count=%d\n", m, n_m, triangles, Ec, Ic);

            EflattenX = malloc(sizeof(double) * Ec);
            EflattenY = malloc(sizeof(double) * Ec);
            IflattenX = malloc(sizeof(double) * Ic);
            IflattenY = malloc(sizeof(double) * Ic);

            //avoid definition of a complex struct
            int i;
            int ind = 0;
            for(i = 0; i < Ec; i++){
                EflattenX[ind] = E[i].x;
                EflattenY[ind] = E[i].y;
                ind++;
                //printf("(%f, %f)\n", E[i].x, E[i].y);
            }
		    compute_counts_and_displs(&sendcounts[0], &displs[0], processes, Ec);
            //int jj;
            //for(jj = 0; jj < processes; jj++){
                //printf("[j: %d, #: %d, from: %d ]", jj, sendcounts[jj], displs[jj]);
            //}
        }

        int maxE = (n_m >= 1) ? (3 * n_m - 3) : 0;
        //here must happen the scatter magic
        int local_length = compute_local_length(processes, rank, maxE);
        Point* local_E = malloc(sizeof(Point) * local_length);
        /*double* local_EX = malloc(sizeof(double) * local_length);
        double* local_EY = malloc(sizeof(double) * local_length);
        MPI_Scatterv(EflattenX, sendcounts, displs, MPI_DOUBLE, local_EX, local_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(EflattenY, sendcounts, displs, MPI_DOUBLE, local_EY, local_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);*/
        MPI_Scatterv(E, sendcounts, displs, MPI_POINT, local_E, local_length, MPI_POINT, 0, MPI_COMM_WORLD);

        double local_sum = 0;
        int indx;
        for(indx = 0; indx < local_length; indx++){
            local_sum += f(local_E[indx].x, local_E[indx].y);
            //printf("[%d-%f-%f]\n", rank, local_EX[indx], local_EY[indx]);
        }
        MPI_Reduce(&local_sum, &sumE, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //free(local_EX);
        //free(local_EY);
        free(local_E);

        if(rank == 0){
            int i;
            int ind = 0;
            for(i = 0; i < Ic; i++){
                IflattenX[ind] = I[i].x;
                IflattenY[ind] = I[i].y;
                ind++;
                //printf("(%f, %f)\n", I[i].x, I[i].y);
            }
		    compute_counts_and_displs(&sendcounts[0], &displs[0], processes, Ic);
        }

        int total_pts = (n_m + 1) * (n_m + 2) / 2;
        int maxI = total_pts - 3 * n_m;
        if(maxI < 0){
            maxI = 0;
        }
        local_length = compute_local_length(processes, rank, maxI);
        Point* local_I = malloc(sizeof(Point) * local_length);
        //double* local_IX = malloc(sizeof(double) * local_length);
        //double* local_IY = malloc(sizeof(double) * local_length);
        //MPI_Scatterv(IflattenX, sendcounts, displs, MPI_DOUBLE, local_IX, local_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Scatterv(IflattenY, sendcounts, displs, MPI_DOUBLE, local_IY, local_length, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(I, sendcounts, displs, MPI_POINT, local_I, local_length, MPI_POINT, 0, MPI_COMM_WORLD);

        local_sum = 0;
        for(indx = 0; indx < local_length; indx++){
            local_sum += f(local_I[indx].x, local_I[indx].y);
            //printf("[%d-%f-%f]\n", rank, local_EX[indx], local_EY[indx]);
        }
        MPI_Reduce(&local_sum, &sumI, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //free(local_IX);
        //free(local_IY);
        free(local_I);
    
        if(rank == 0){
            sumE*=3.0;
            sumI*=6.0;
            //last formula computation
            acc = (area / (3.0 * triangles)) * (sumC + sumE + sumI);
            R[m][0] = acc;

            printf("R[%d][0] = %f\n", m, R[m][0]);
            free(E);
            free(I);
            free(EflattenX);
            free(EflattenY);
            printf("\n");       
        }

    }

    if(rank == 0){
        //Richardson Extraploation
        printf("R[0]: %f\n", R[0][0]);
        int m;
        for(m = 1; m < MaxLevel; m++){
            //compute R[m][1], R[m][2], ..., R[m][m]
            printf("R[%d]: %f ", m, R[m][0]);
            int k;
            for(k = 1; k <= m; k++){    
                R[m][k] = R[m][k-1] + (-R[m-1][k-1] + R[m][k-1]) / (pow(2, k) - 1); //on the paper sign into the () are opposite
                printf("%f ", R[m][k]);
            }
            printf("\n");
        }

        //deallocation of R
        size_t i;
        for(i = 0; i < MaxLevel; i++){
            free(R[i]);
        }
        free(R);

        double endTime = MPI_Wtime();
        printf("Time elapsed %f, %d processors, complexity %d", endTime - startTime, processes, ITER);
    }

	MPI_Type_free(&MPI_POINT);
    MPI_Finalize();



    return 0;
}