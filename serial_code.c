#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define ITER 100

double f(double x, double y) {

    double result = 0;

    for(int i=0; i < ITER; i++){
        result += pow(x, 2) + pow(y, 2) + 1;
    }

    return 1 / (result);
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
        for (int t = 1; t <= nm - 1; ++t) { //compute all points on each edge (0...nm) [I have nm+1 points][0 is A so I exclude it][nm is B so I exclude it]
            double u = (double)t / nm; //contribution of each point
            
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
    for (int ia = 1; ia <= nm - 2; ++ia) {
        /* ib can range from 1 to nm-1-ia, such that ic = nm-ia-ib >=1 */
        for (int ib = 1; ib <= nm - 1 - ia; ++ib) {
            
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



int main(int argc, char* argv[]){
    
    //paremeters chosen
    int MaxLevel = 8;

    Point A = {2.0, 1.0};
    Point B = {4.0, 2.5};
    Point C = {6.0, 0.75};


    //INPUT PREPARATION
    double** R = malloc(MaxLevel * sizeof(double*));
    for(size_t i = 0; i < MaxLevel; i++){
        R[i] = malloc(MaxLevel * sizeof(double));
    }

    double area = 0.5 * fabs(A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
    printf("Area is %f\n", area);

    for(int m = 0; m < MaxLevel; m++){
        //compute R[0][0], R[1][0], ..., R[m][0], ..., R[MaxLevel][0]
        double acc = 0;
        double n_m = (double) (1 << m);
        int triangles = (int) (pow(n_m, 2));

        double sumC = 0, sumE = 0, sumI = 0; //C, Em, Im

        //first sum member
        sumC += f(A.x, A.y) + f(B.x, B.y) + f(C.x, C.y);

        //edge and interior nodes computation

        Point *E = NULL, *I = NULL;
        int Ec = 0, Ic = 0;
        generate_Em_Im(A, B, C, (int) n_m, &E, &Ec, &I, &Ic);
        printf("m=%d  nm=%f triangles=%d edge_count (excluding vertices)=%d  interior_count=%d\n", m, n_m, triangles, Ec, Ic);

        //second sum member
        for (int i = 0; i < Ec; ++i){
            sumE += f(E[i].x, E[i].y);
        }
        sumE *= 3;

        for (int i = 0; i < Ic; ++i){
            sumI += f(I[i].x, I[i].y);
        }
        sumI*= 6;

        //last formula computation
        acc = (area / (3.0 * triangles)) * (sumC + sumE + sumI);
        R[m][0] = acc;

        printf("R[%d][0] = %f\n", m, R[m][0]);
        free(E);
        free(I);
        printf("\n");
    }

    //Richardson Extraploation
    printf("R[0]: %f\n", R[0][0]);
    for(int m = 1; m < MaxLevel; m++){
        //compute R[m][1], R[m][2], ..., R[m][m]
        printf("R[%d]: %f ", m, R[m][0]);
        for(int k = 1; k <= m; k++){
            R[m][k] = R[m][k-1] + (-R[m-1][k-1] + R[m][k-1]) / (pow(2, k) - 1); //on the paper sign into the () are opposite
            printf("%f ", R[m][k]);
        }
        printf("\n");
    }

    //deallocation of R
    for(size_t i = 0; i < MaxLevel; i++){
        free(R[i]);
    }
    free(R);

    return 0;
}

