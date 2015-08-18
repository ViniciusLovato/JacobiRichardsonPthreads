#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void MAandMBdividebydiagonal(int J_ORDER, double** MA, double MB[]);
double* MAstartimesXk(int J_ORDER, double** MA, double Xk[]);
double* sub_MBMAstartimesXk(int J_ORDER, double MB[], double MAstartimes[]);
void JacobiRichardson(int J_ORDER, int J_ITE_MAX, double* J_ERROR, double** MA, double MB[]);
double get_error(int J_ORDER, double X[], double Xk[]);
double max_arrayval(int J_ORDER, double ar[]);
void printArray(int J_ORDER, double ar[]);
void printMatrix(int J_ORDER, double** Matrix);
void readStuffFromFile(FILE* pFile, int* J_ORDER, int* J_ROW_TEST, double* J_ERROR, int* J_ITE_MAX, double*** MA, double** MB);
void cleanup(FILE* pFile, int J_ORDER, double** MA, double MB[]);

int main(int argc, char* argv[]) { 
	FILE* pFile;
    pFile = fopen (argv[1] , "r");
    int J_ORDER, J_ROW_TEST, J_ITE_MAX;
    double J_ERROR;
    double** MA = NULL;
    double* MB = NULL;
    readStuffFromFile(pFile, &J_ORDER, &J_ROW_TEST, &J_ERROR, &J_ITE_MAX, &MA, &MB);
    MAandMBdividebydiagonal(J_ORDER, MA, MB);
    JacobiRichardson(J_ORDER, J_ITE_MAX, &J_ERROR, MA, MB);
    cleanup(pFile, J_ORDER, MA, MB);
}

void MAandMBdividebydiagonal(int J_ORDER, double** MA, double MB[]) {
	
    for (int i = 0; i < J_ORDER; ++i) {
        for (int j = 0; j < J_ORDER; ++j) {
            if (i != j) {
                MA[i][j] /= MA[i][i];//divide todo mundo da linha pela diagonal
            }
        }
        MB[i] /= MA[i][i];//divide a matriz B pela diagonal
        MA[i][i] = 0;//faz a diagonal ser 0
    }
}

void JacobiRichardson(int J_ORDER, int J_ITE_MAX, double* J_ERROR, double** MA, double MB[]) {
    int iterations = 0;
    double e = 100;
    double* X;
    double* MAstartimes;
    double* Xk = (double*) malloc(sizeof(double) * J_ORDER);
    memcpy(Xk, MB, sizeof(double) * J_ORDER);

    while (e > *J_ERROR && iterations < J_ITE_MAX) {
        ++iterations;
        MAstartimes = MAstartimesXk(J_ORDER, MA, Xk);
        X = sub_MBMAstartimesXk(J_ORDER, MB, MAstartimes);
        e = get_error(J_ORDER, X, Xk);
        memcpy(Xk, X, sizeof(double) * J_ORDER);

    }
	printf("iterations: %d\n", iterations);
	printArray(J_ORDER, X);
    free(X);
    free(Xk);
    free(MAstartimes);
}

double* MAstartimesXk(int J_ORDER, double** MA, double Xk[]) {
    double* vect = (double*) calloc(sizeof(double), J_ORDER);
    for (int i = 0; i < J_ORDER; ++i) {
        for (int j = 0;  j < J_ORDER; ++j) {
            vect[i] += MA[i][j] * Xk[j];
        }
    }
    return vect;
}

double* sub_MBMAstartimesXk(int J_ORDER, double MB[], double MAstartimes[]) {
    double* vect = (double*) malloc(sizeof(double) * J_ORDER);
    for (int i = 0; i < J_ORDER; ++i) {
        vect[i] = MB[i] - MAstartimes[i];
    }
    return vect;

}

double get_error(int J_ORDER, double X[], double Xk[]) {
    double* dif = (double*) malloc(sizeof(double) * J_ORDER);
    for (int i = 0; i < J_ORDER; ++i) {
        dif[i] = fabs(X[i]-Xk[i]);
    }

    double err = max_arrayval(J_ORDER, dif)/max_arrayval(J_ORDER, X);
    free(dif);
    return err;

}

double max_arrayval(int J_ORDER, double ar[]) {
    double maxValue = ar[0];

    for (int i = 1; i < J_ORDER; ++i) {
        if (ar[i] > maxValue) {
            maxValue = ar[i];
        }
    }

    return maxValue;
}

void printArray(int J_ORDER, double ar[]) {
    for (int i = 0; i < J_ORDER; ++i) {
        printf("%f ", ar[i]);
    }
    printf("\n");
}

void printMatrix(int J_ORDER, double** Matrix) {
    for (int i = 0; i < J_ORDER; ++i) {
        for (int j = 0; j < J_ORDER; ++j) {
            printf("%f ", Matrix[i][j]);
        }
        printf("\n");
    }
}

void cleanup(FILE* pFile, int J_ORDER, double** MA, double* MB) {
    for (int i = 0; i < J_ORDER; ++i){
        free(MA[i]);
    }
    free(MA);
    free(MB);
    fclose(pFile);
}

void readStuffFromFile(FILE* pFile, int* J_ORDER, int* J_ROW_TEST, double* J_ERROR, int* J_ITE_MAX, double*** MA, double** MB) {
	fscanf(pFile, "%d%d%lf%d", J_ORDER,J_ROW_TEST,J_ERROR,J_ITE_MAX);

    *MA = (double**) malloc(sizeof(double*) * (*J_ORDER));
    for (int i = 0; i < (*J_ORDER); ++i){
        MA[0][i] = (double*) malloc(sizeof(double) * (*J_ORDER));
    }
    for (int i = 0; i < (*J_ORDER); ++i) {
        for (int j = 0; j < (*J_ORDER); ++j) {
            fscanf(pFile, "%lf", &MA[0][i][j]);
        }
    }
    *MB = (double*) malloc(sizeof(double) * (*J_ORDER));
    for (int i = 0; i < (*J_ORDER); ++i) {
       fscanf(pFile, "%lf", &MB[0][i]);
    }
}
