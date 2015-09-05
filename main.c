#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/*modo de uso: programa nome_do_arquivo_de_entrada
Eg: main teste.txt

pra compilar:
gcc main.c -o main -std=c99
*/


/*divide matriz A e matriz B(aka vetor B) pela diagonal da Matriz A, e faz diagonal de Matriz A ser 0. Isso faz a matriz A virar matriz A*(matriz A modificada, aka MAstar) e a matriz B virar matriz B*(matriz B modificada, aka MBstar). Ver pdf explicando o metodo de Jacobi. Recebe como parametro a ordem da matriz e as matrizes*/
void MAandMBdividebydiagonal(int J_ORDER, double** MA, double MB[], double** MAstar, double MBstar[]);
/*multiplica a matriz A* por Xk. Ver pdf explicando o metodo de Jacobi.Recebe como parametro a ordem da matriz e as matrizes e retorna a matriz(vetor) resultante*/
void MAstartimesXk(int J_ORDER, double** MAstar, double MAstartimes[], double Xk[]);
/*subtrai matriz B* por matrizA* multiplicado por Xk. Faz parte da formula da iteratividade do metodo de Jacobi. Recebe como parametro a ordem da matriz, a matriz MB*, a matriz A* multiplicado por Xk(MAstartimes) e o vetor Xk+1*/
void sub_MBMAstartimesXk(int J_ORDER, double MBstar[], double MAstartimes[], double X[]);
/*aplica o metodo de Jacobi propriamente dito. Recebe a ordem da matriz, numero maximo de iteracoes permitido, linha pra testar, erro, matriz MA* e matriz MB* */
void JacobiRichardson(int J_ORDER, int J_ITE_MAX, int J_ROW_TEST, double* J_ERROR, double** MA, double MB[], double** MAstar, double MBstar[]);
/*calcula erro pela formula do pdf. Recebe ordem da matriz, o vetor de X atual(X) e o vetor de X anterior(Xk)*/
double get_error(int J_ORDER, double X[], double Xk[]);
/*funcao auxiliar para retornar o maior valor de um vetor, recebe ordem do vetor e o vetor. Retorna o maior valor encontrado. Eh usado para calcular o erro*/
double max_arrayval(int J_ORDER, double ar[]);
/*imprime vetor. Recebe ordem do vetor e o vetor. Funcao para debug*/
void printArray(int J_ORDER, double ar[]);
/*imprime matriz. Recebe ordem da matriz e a matriz. Funcao para debug*/
void printMatrix(int J_ORDER, double** Matrix);
/*pega as informacoes do arquivo de acordo com a especificacao do trabalho e salva nos parametros passado como argumento da funcao.
Os parametros tem exatamente o mesmo nome da descricao do trabalho, todos devem ser passados por referencia, inclusive as matrizes, porque elas sao mallocadas dentro da funcao*/
void readStuffFromFile(FILE* pFile, int* J_ORDER, int* J_ROW_TEST, double* J_ERROR, int* J_ITE_MAX, double*** MA, double** MB);
/*para dar free nas matrizes MA, MB,  MA* e MB* e fechar arquivo de entrada*/
void cleanup(FILE* pFile, int J_ORDER, double** MA, double* MB, double** MAstar, double MBstar[]);
/*para calcular resultado em determinada linha. Recebe vetor X final, matriz MA original e a linha que deve ser usada*/
double result_row(double* X, double** MA, int J_ROW_TEST, int J_ORDER);

int main(int argc, char* argv[]) { 
	FILE* pFile;//ponteiro pro arquivo de entrada
    pFile = fopen (argv[1] , "r");//abre arquivo de entrada para leitura somente
    int J_ORDER, J_ROW_TEST, J_ITE_MAX;//variaveis que vao armazenar valores do arquivo. Ordem da matriz, linha de teste, maximo de iteracao
    double J_ERROR;//armazena erro
    double** MA = NULL;//matriz A original
    double* MB = NULL;//matriz(vetor) B original
    readStuffFromFile(pFile, &J_ORDER, &J_ROW_TEST, &J_ERROR, &J_ITE_MAX, &MA, &MB);//le o arquivo e armazena nas variaveis os valores
	double** MAstar = (double**) malloc(sizeof(double*) * J_ORDER);//matriz A*(matriz A modificada), gerada na funcao MAandMBdividebydiagonal
	for (int i = 0; i < J_ORDER; ++i) MAstar[i] = (double*) malloc(sizeof(double) * J_ORDER);
	double* MBstar = (double*) malloc(sizeof(double) * J_ORDER);//matriz B*(matriz B modificada), gerada na funcao MAandMBdividebydiagonal
    MAandMBdividebydiagonal(J_ORDER, MA, MB, MAstar, MBstar);//divide MA e MB pelos valores da diagonal de MA, produzindo MA* e MB*
    JacobiRichardson(J_ORDER, J_ITE_MAX, J_ROW_TEST, &J_ERROR, MA, MB, MAstar, MBstar);//aplica metodo de Jacobi
    cleanup(pFile, J_ORDER, MA, MB, MAstar, MBstar);//da free nos malloc de MA e MB, e fecha arquivo de entrada
}

void MAandMBdividebydiagonal(int J_ORDER, double** MA, double MB[], double** MAstar, double MBstar[]) {
    for (int i = 0; i < J_ORDER; ++i) {
        for (int j = 0; j < J_ORDER; ++j) {
				MAstar[i][j] = MA[i][j]/MA[i][i];//divide todo mundo da linha pela diagonal
        }
		MBstar[i] = MB[i]/MA[i][i];//divide a matriz B pela diagonal
        MAstar[i][i] = 0;//faz a diagonal ser 0
    }
}

void JacobiRichardson(int J_ORDER, int J_ITE_MAX, int J_ROW_TEST, double* J_ERROR, double** MA, double MB[], double** MAstar, double MBstar[]) {    
	int iterations = 1;
    double e = 100;
    double* X = (double*) malloc(sizeof(double) * J_ORDER);//X eh o Xk+1
    double* MAstartimes = (double*) malloc(sizeof(double) * J_ORDER);
    double* Xk = (double*) malloc(sizeof(double) * J_ORDER);//aloca Xk para armazenar sempre a iteracao anterior.
    memcpy(Xk, MBstar, sizeof(double) * J_ORDER);//inicialmente Xk(na primeira iteracao, o X0) recebe os valores do vetor B*

    while (e > *J_ERROR && iterations < J_ITE_MAX) {//enquanto for maior que erro permitido e nao atingir numero maximo de iteracao
        ++iterations;//incremente numero de iteracoes
        MAstartimesXk(J_ORDER, MAstar, MAstartimes, Xk);//calcula MA* vezes Xk     
		sub_MBMAstartimesXk(J_ORDER, MBstar, MAstartimes, X);//subtrai MB* do resultado anterior. X eh o Xk+1
        e = get_error(J_ORDER, X, Xk);//calcula erro
        memcpy(Xk, X, sizeof(double) * J_ORDER);//atualiza Xk

    }

	//imprime numero de iteracoes e resultado em determinada linha
	printf("iterations: %d\nRowTest: %d => [%f] =? %f\n", iterations, J_ROW_TEST, result_row(X, MA, J_ROW_TEST, J_ORDER), MB[J_ROW_TEST]);
    free(X);
    free(Xk);
    free(MAstartimes);
}

double result_row(double* X, double** MA, int J_ROW_TEST, int J_ORDER) {
	double result = 0.0;
	for (int i = 0; i < J_ORDER; ++i) {
		result += MA[J_ROW_TEST][i] * X[i];
	}

	return result;
}

void MAstartimesXk(int J_ORDER, double** MAstar, double MAstartimes[], double Xk[]) {
	//multiplicacao de matriz por vetor
    for (int i = 0; i < J_ORDER; ++i) {
		MAstartimes[i] = 0.0;//zera vetor primeiro
        for (int j = 0;  j < J_ORDER; ++j) {
            MAstartimes[i] += MAstar[i][j] * Xk[j];
        }
    }
}

void sub_MBMAstartimesXk(int J_ORDER, double MBstar[], double MAstartimes[], double X[]) {
	//subtrai MB* de MA* vezes Xk e armazena em Xk+1
    for (int i = 0; i < J_ORDER; ++i) {
        X[i] = MBstar[i] - MAstartimes[i];
    }

}

double get_error(int J_ORDER, double X[], double Xk[]) {
	//calcula erro seguindo formula do pdf
    double* dif = (double*) malloc(sizeof(double) * J_ORDER);
    for (int i = 0; i < J_ORDER; ++i) {
        dif[i] = fabs(X[i]-Xk[i]);
    }

    double err = max_arrayval(J_ORDER, dif)/max_arrayval(J_ORDER, X);
    free(dif);
    return err;

}

double max_arrayval(int J_ORDER, double ar[]) {
	//retorna valor maximo de vetor
    double maxValue = fabs(ar[0]);
	double nextabsolute;

    for (int i = 1; i < J_ORDER; ++i) {
		nextabsolute = fabs(ar[i]);
        if (nextabsolute > maxValue) {
            maxValue = nextabsolute;
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

void cleanup(FILE* pFile, int J_ORDER, double** MA, double* MB, double** MAstar, double MBstar[]) {
    for (int i = 0; i < J_ORDER; ++i){
        free(MA[i]);
		free(MAstar[i]);
    }
    free(MA);
	free(MAstar);
    free(MB);
	free(MBstar);
    fclose(pFile);
}

void readStuffFromFile(FILE* pFile, int* J_ORDER, int* J_ROW_TEST, double* J_ERROR, int* J_ITE_MAX, double*** MA, double** MB) {
	fscanf(pFile, "%d%d%lf%d", J_ORDER,J_ROW_TEST,J_ERROR,J_ITE_MAX);
	//matriz MA eh passado por referencia pra poder dar malloc dentro da funcao. Eh feio, mas eh assim em C.
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
