#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

/*modo de uso: programa nome_do_arquivo_de_entrada numero_threads
Eg: parallel teste.txt 2

pra compilar:
gcc parallel.c -o parallel -std=c99 -pthread
*/

/*cria as threads. Foi necessario usar macro porque sem ela da erro de "dereference void* pointer"*/
#define create_threads(start_routine, arg, thr, nthreads); \
for (int i = 0, rc = 0; i < nthreads; ++i) {\
	if ((rc = pthread_create(&thr[i], NULL, start_routine, &arg[i]))) {\
		fprintf(stderr, "error: pthread_create, rc: %d\n", rc);\
		return EXIT_FAILURE;\
	}\
}
		
typedef struct {
	int J_ORDER;
	int begin;
	int end;
	double** MA;
	double* MB;
	double** MAstar;
 	double* MBstar;
} MAandMBdividebydiagonal_args;

typedef struct {
	int J_ORDER;
	int begin;
	int end;
	double** MAstar;
	double* MAstartimes;
	double* Xk;
} MAstartimesXk_args;


/*espera as threads terminarem antes de continar a main*/
void wait_threads(pthread_t thr[], int nthreads);
/*seta os vetores de comeco e fim do trabalho de cada thread. Cada thread vai trabalhar em um bloco da matriz*/
void set_begin_end(int begin[], int end[], int nthreads, int total_work_for_each_thread, int total_work_for_first_thread);

/*seta os argumentos da funcao MAandMBdividebydiagonal. Recebe como parametro o vetor de argumentos para setar(arguments) e as variaveis que vao ser usadas para isso(J_ORDER, begin, etc)*/
void setMAandMBdividebydiagonal_args(MAandMBdividebydiagonal_args arguments[], int J_ORDER, int begin[], int end[], double** MA, double MB[], double** MAstar, double MBstar[], int nthreads);
/*divide matriz A e matriz B(aka vetor B) pela diagonal da Matriz A, e faz diagonal de Matriz A ser 0. Isso faz a matriz A virar matriz A*(matriz A modificada, aka MAstar) e a matriz B virar matriz B*(matriz B modificada, aka MBstar). Ver pdf explicando o metodo de Jacobi. Recebe como parametro a ordem da matriz e as matrizes*/
void* MAandMBdividebydiagonal(void* args);

/*seta os argumentos da funcao MAstartimesXk. Recebe como parametro o vetor de argumentos para setar(arguments) e as variaveis que vao ser usadas para isso(J_ORDER, begin, etc)*/
void setMAstartimesXk_args(MAstartimesXk_args arguments[], int J_ORDER, int begin[], int end[], double** MAstar, double* MAstartimes, double Xk[], int nthreads);
/*multiplica a matriz A* por Xk*/
void* MAstartimesXk(void* args);


/*subtrai matriz B* por matrizA* multiplicado por Xk. Faz parte da formula da iteratividade do metodo de Jacobi. Retorna a matriz(vetor) resultante. Recebe como parametro a ordem da matriz, a matriz MB* e a matriz A* multiplicado por Xk(MAstartimes)*/
double* sub_MBMAstartimesXk(int J_ORDER, double MBstar[], double MAstartimes[]);
/*aplica o metodo de Jacobi propriamente dito. Recebe a ordem da matriz, numero maximo de iteracoes permitido, linha pra testar, erro, matriz MA*, matriz MB*, numero de threads, vetor de comeco, fim, e as threads*/
int JacobiRichardson(int J_ORDER, int J_ITE_MAX, int J_ROW_TEST, double* J_ERROR, double** MA, double MB[], double** MAstar, double MBstar[], int nthreads, int begin[], int end[], pthread_t thr[]);
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

	int err;//pra caso nao aloque as threads
	FILE* pFile;//ponteiro pro arquivo de entrada
    pFile = fopen (argv[1] , "r");//abre arquivo de entrada para leitura somente
    int nthreads = atoi(argv[2]);//armazena numero de threads(nao conta com a thread main)

	int J_ORDER, J_ROW_TEST, J_ITE_MAX;//variaveis que vao armazenar valores do arquivo. Ordem da matriz, linha de teste, maximo de iteracao
    double J_ERROR;//armazena erro
    double** MA = NULL;//matriz A original
    double* MB = NULL;//matriz(vetor) B original
    readStuffFromFile(pFile, &J_ORDER, &J_ROW_TEST, &J_ERROR, &J_ITE_MAX, &MA, &MB);//le o arquivo e armazena nas variaveis os valores

	int total_work_for_each_thread = J_ORDER/nthreads;//armazena quantidade de linhas de uma matriz que cada thread vai cuidar
	int total_work_for_first_thread = total_work_for_each_thread + J_ORDER%nthreads;//nem sempre todas as threads vao ter uma quantidade igual de trabalho. Nesse caso, a primeira thread fica com um pouco mais de trabalho pra fazer
	int begin[nthreads], end[nthreads];//cria vetor de onde comeca e termina trabalho de cada thread. Por exemplo, se a matriz for 4x4 e tiver 2 threads, a primeira thread trabalha nas linhas 0 e 1, a segunda trabalha nas linhas 2 e 3. Cada thread fica com um bloco de tamanho igual da matriz para trabalhar. Se nao der pra dividir igualmente, a primeira thread pega um pouco +.
	pthread_t thr[nthreads];//cria vetor de threads
	set_begin_end(begin, end, nthreads, total_work_for_each_thread, total_work_for_first_thread);//seta os vetores de comeco e fim do trabalho de cada thread

	double** MAstar = (double**) malloc(sizeof(double*) * J_ORDER);//matriz A*(matriz A modificada), gerada na funcao MAandMBdividebydiagonal
	for (int i = 0; i < J_ORDER; ++i) MAstar[i] = (double*) malloc(sizeof(double) * J_ORDER);
	double* MBstar = (double*) malloc(sizeof(double) * J_ORDER);//matriz B*(matriz B modificada), gerada na funcao MAandMBdividebydiagonal


	MAandMBdividebydiagonal_args MAandMBdividebydiagonal_arguments[nthreads];//cria vetor de argumentos pra funcao MAandMBdividebydiagonal
	setMAandMBdividebydiagonal_args(MAandMBdividebydiagonal_arguments, J_ORDER, begin, end, MA, MB, MAstar, MBstar, nthreads);//seta os argumentos da funcao MAandMBdividebydiagonal
	create_threads(MAandMBdividebydiagonal, MAandMBdividebydiagonal_arguments, thr, nthreads);//cria as threads
	wait_threads(thr, nthreads);//espera as threads acabarem antes de continuar
		
   
    err = JacobiRichardson(J_ORDER, J_ITE_MAX, J_ROW_TEST, &J_ERROR, MA, MB, MAstar, MBstar, nthreads, begin, end, thr);//aplica metodo de Jacobi
    cleanup(pFile, J_ORDER, MA, MB, MAstar, MBstar);//da free nos malloc de MA e MB, e fecha arquivo de entrada
	return err;
}

void wait_threads(pthread_t thr[], int nthreads) {
	for (int i = 0; i < nthreads; ++i) {
		pthread_join(thr[i], NULL);
	}
}

void set_begin_end(int begin[], int end[], int nthreads, int total_work_for_each_thread, int total_work_for_first_thread) {
	begin[0] = 0;
	end[0]   = total_work_for_first_thread - 1;
	for (int i = 1; i < nthreads; ++i) {
		begin[i] = end[i - 1] + 1;
		end[i] = begin[i] + total_work_for_each_thread - 1;
	}
}

void setMAandMBdividebydiagonal_args(MAandMBdividebydiagonal_args arguments[], int J_ORDER, int begin[], int end[], double** MA, double MB[], double** MAstar, double MBstar[], int nthreads) {
	for (int i = 0; i < nthreads; ++i) {
		arguments[i].J_ORDER = J_ORDER;
		arguments[i].begin   = begin[i];
		arguments[i].end     = end[i];
		arguments[i].MA      = MA;
		arguments[i].MB      = MB;
		arguments[i].MAstar  = MAstar;
		arguments[i].MBstar  = MBstar;
	}
}

void* MAandMBdividebydiagonal(void* args) {
	MAandMBdividebydiagonal_args* arguments = args;
	int J_ORDER 	= arguments->J_ORDER;
	int begin   	= arguments->begin;
	int end     	= arguments->end;
	double** MA 	= arguments->MA;
	double* MB 		= arguments->MB;
	double** MAstar = arguments->MAstar;
 	double* MBstar 	= arguments->MBstar;

	for (int i = begin; i <= end; ++i) {
        for (int j = 0; j < J_ORDER; ++j) {
				MAstar[i][j] = MA[i][j]/MA[i][i];//divide todo mundo da linha pela diagonal
        }
		MBstar[i] = MB[i]/MA[i][i];//divide a matriz B pela diagonal
        MAstar[i][i] = 0;//faz a diagonal ser 0
    }

	pthread_exit(NULL);
}

void setMAstartimesXk_args(MAstartimesXk_args arguments[], int J_ORDER, int begin[], int end[], double** MAstar, double* MAstartimes, double Xk[], int nthreads) {
	for (int i = 0; i < nthreads; ++i) {
		arguments[i].J_ORDER 	  = J_ORDER;
		arguments[i].begin   	  = begin[i];
		arguments[i].end     	  = end[i];
		arguments[i].MAstar  	  = MAstar;
		arguments[i].MAstartimes  = MAstartimes;
		arguments[i].Xk           = Xk;
	}
}

void* MAstartimesXk(void* args) {
	MAstartimesXk_args* arguments = args;
	int J_ORDER 		= arguments->J_ORDER;
	int begin   		= arguments->begin;
	int end     		= arguments->end;
	double** MAstar 	= arguments->MAstar;
	double* MAstartimes = arguments->MAstartimes;
 	double* Xk 			= arguments->Xk;

	//multiplicacao de matriz por vetor
    for (int i = begin; i <= end; ++i) {
		MAstartimes[i] = 0.0;//zera vetor primeiro
        for (int j = 0;  j < J_ORDER; ++j) {
            MAstartimes[i] += MAstar[i][j] * Xk[j];
        }
    }
	
    pthread_exit(NULL);
}


int JacobiRichardson(int J_ORDER, int J_ITE_MAX, int J_ROW_TEST, double* J_ERROR, double** MA, double MB[], double** MAstar, double MBstar[], int nthreads, int begin[], int end[], pthread_t thr[]) {    	
	int iterations = 1;
    double e = 100;
    double* X;//X eh o Xk+1
    double* MAstartimes = (double*) malloc(sizeof(double) * J_ORDER);
    double* Xk = (double*) malloc(sizeof(double) * J_ORDER);//aloca Xk para armazenar sempre a iteracao anterior.

    memcpy(Xk, MBstar, sizeof(double) * J_ORDER);//inicialmente Xk(na primeira iteracao, o X0) recebe os valores do vetor B*
	MAstartimesXk_args MAstartimesXk_arguments[nthreads];//cria vetor de argumentos pra funcao MAstartimesXk
	
	setMAstartimesXk_args(MAstartimesXk_arguments, J_ORDER, begin, end, MAstar, MAstartimes, Xk, nthreads);//seta os argumentos da funcao MAstartimesXk
    while (e > *J_ERROR && iterations < J_ITE_MAX) {//enquanto for maior que erro permitido e nao atingir numero maximo de iteracao
        ++iterations;//incremente numero de iteracoes
		create_threads(MAstartimesXk, MAstartimesXk_arguments, thr, nthreads);//cria as threads
		wait_threads(thr, nthreads);//espera as threads acabarem antes de continuar
		X = sub_MBMAstartimesXk(J_ORDER, MBstar, MAstartimes);//subtrai MB* do resultado anterior. X eh o Xk+1
        e = get_error(J_ORDER, X, Xk);//calcula erro
        memcpy(Xk, X, sizeof(double) * J_ORDER);//atualiza Xk
    }
	
	//imprime numero de iteracoes e resultado em determinada linha
	printf("iterations: %d\nRowTest: %d => [%f] =? %f\n", iterations, J_ROW_TEST, result_row(X, MA, J_ROW_TEST, J_ORDER), MB[J_ROW_TEST]);
    free(X);
    free(Xk);
    free(MAstartimes);
	return 0;
}

double result_row(double* X, double** MA, int J_ROW_TEST, int J_ORDER) {
	double result = 0.0;
	for (int i = 0; i < J_ORDER; ++i) {
		result += MA[J_ROW_TEST][i] * X[i];
	}

	return result;
}

double* sub_MBMAstartimesXk(int J_ORDER, double MBstar[], double MAstartimes[]) {
	//subtrai MB* de MA* vezes Xk
    double* vect = (double*) malloc(sizeof(double) * J_ORDER);
    for (int i = 0; i < J_ORDER; ++i) {
        vect[i] = MBstar[i] - MAstartimes[i];
    }
    return vect;

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
