#include <stdio.h>
#include <stdlib.h>


double** readA(int);
double* readB(int);
void printMatrix(double**, int);
void freeAll(double**, int,  double*);

/*
 * This algorithim calculates Ax = B
 * using Jacobi Richardson method
 *
 */
int main(){
  
  	// matrix order
	int J_ORDER = 0;

	// matrix row to test
	int J_ROW_TEST = 0;	

	// error tolerance
	double J_ERROR = 0;

	// maximum number of iterations allowed
	int J_ITE_MAX = 0;

	// Matrix A
	double **A;

	// Array B
	double *b;

	scanf("%d", &J_ORDER);
	scanf("%d", &J_ROW_TEST);
	scanf("%lf", &J_ERROR);
	scanf("%d", &J_ITE_MAX);

    A = readA(J_ORDER);
    b = readB(J_ROW_TEST);

    //printMatrix(A, J_ORDER);

    freeAll(A, J_ORDER,  b);
    return 0;
}

/*
 * Reads the matrix A. Ax = B
 * @return Matrix A
 *
 */
double** readA(int n){
    
    double** A = NULL;
    int i, j = 0;
    // 
    A = (double**) malloc(sizeof(double*) * n); 
    for (i = 0; i < n; i++){
        A[i] = (double*) malloc(sizeof(double) * n);
    }

    // 
    for (i = 0; i < n; i++){
		for (j = 0; j < n; j++){
			scanf("%lf", &A[i][j]);
		}
	}
	
	return A;
}

/*
 * Reads the array b. Ax = b
 * @return array b
 *
 */
double * readB(int n){
	
	double *b = NULL;
    int i = 0;

	b = (double*) malloc(sizeof(double) * n);

	for(i = 0; i < n; i++){
		scanf("%lf", &b[i]);
	}

	return b;
}

/*
 * prints matrix
 * @return 
 *
 */
void printMatrix(double** M, int n){
    
    int i, j = 0;

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            printf("%lf ", M[i][j]);
        }
        printf("\n");
    }
    return;
}

/*
 * free all allocated memory used by matrix A and array B 
 * @param A 
 * @param b
 *
 */
void freeAll(double** A, int n,  double* b){

    int i = 0;
    free(b);
   
    for(i = 0; i < n; i++){
        free(A[i]);
    }

    free(A);
   
    return;
}
