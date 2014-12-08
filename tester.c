#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <papi.h>

#define NUM_EVENTS 2
#define SIZE 10


void createValidMatrices(float a[SIZE][SIZE], float b[SIZE][SIZE], float res[SIZE][SIZE]){
        int i,j;
        srand((unsigned)time(NULL));
        for(i = 0; i < SIZE; i++){
                for(j = 0; j < SIZE; j++){
                        a[i][j] = ((float)rand()) / ((float)RAND_MAX + 1);
                        b[i][j] = 1;
                        res[i][j] = 0;
                }
        }
}

void createMatrices(float a[SIZE][SIZE], float b[SIZE][SIZE],float res[SIZE][SIZE]){
        int i, j;
        srand((unsigned)time(NULL));
        for(i = 0; i < SIZE; i++){
                for(j = 0; j < SIZE; j++){
                        a[i][j] = ((float)rand()) / ((float)RAND_MAX + 1);
                        b[i][j] = ((float)rand()) / ((float)RAND_MAX + 1);
                        res[i][j] = 0;
                }
        }
}

void multMatrices(float a[SIZE][SIZE], float b[SIZE][SIZE], float res[SIZE][SIZE]){
        int i,j,k;
        for(j = 0; j < SIZE; j++){
                for(i = 0; i < SIZE; i++){
                        res[i][j] = 0.0;
                        for(k = 0; k < SIZE; k++){
                                res[i][j] += a[i][k] * b[k][j];
                        }
                }
        }
}

int validateProgram1(float res[SIZE][SIZE]){
        int i,j, result=1;
        float ref=0, aux;
        for(i = 0; i < SIZE; i++){
                aux = 0;
                for(j = 0; j < SIZE; j++){
                        if(i==0){
                                ref += res[j][i];
                        }
                        aux += res[j][i];
                }
                if(aux != ref){
                        result = 0;
                }
        }
        return result;
}
int validateProgram2(float res[SIZE][SIZE]){
        int i,j, result=1;
        float ref=0.0, aux;
        for(i = 0; i < SIZE; i++){
                aux = 0.0;
                for(j = 0; j < SIZE; j++){
                        if(i==0){
                                ref += res[i][j];
                        }
                        aux += res[i][j];
                }
                if(aux != ref){
                        result = 0;
                }
        }
        return result;
}

void printMatrix(float m[SIZE][SIZE], char nome[]){
        int i, j;
        printf("---- %s ----\n",nome);
        for(i = 0; i < SIZE; i++){
                for(j = 0; j < SIZE; j++){
                        printf("%.2f ",m[i][j]);
                }
                printf("\n");
        }
        printf("\n");
}

int main(int argc, char *argv[]){
float a[SIZE][SIZE], b[SIZE][SIZE], res[SIZE][SIZE];

int Events[]={PAPI_FP_OPS,PAPI_L3_TCM};
long long values[NUM_EVENTS];

if(argc!=2){
        printf("Argumentos inválido\n USAGE: prog_name FLAG\n");
        return 1;
}


if(!strcmp(argv[1],"val1")){
        createValidMatrices(a,b,res);
        if(PAPI_start_counters(Events,NUM_EVENTS) != PAPI_OK) printf("Nao tem contadores papi\n");
        if(PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK) printf("Erro ao ler evento\n");
		multMatrices(a,b,res);
        if(PAPI_stop_counters(values,NUM_EVENTS) != PAPI_OK) printf("Erro ao parar evento\n");
        printMatrix(a,"A");
        printMatrix(b,"B");
        printMatrix(res, "A*B");
        int val = validateProgram1(res);
        if(val==1){
                printf("PROGRAM IS VAlIDATED for A*B\n");
        }else{
                printf("PROGRAM HAS ERRORS!!!\n");
        }
}

if(!strcmp(argv[1],"val2")){
        createValidMatrices(a,b,res);
        if(PAPI_start_counters(Events,NUM_EVENTS) != PAPI_OK) printf("Nao tem contadores papi\n");
        if(PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK) printf("Erro ao ler evento\n");
		multMatrices(b,a,res);
        if(PAPI_stop_counters(values,NUM_EVENTS) != PAPI_OK) printf("Erro ao parar evento\n");
        printMatrix(a,"A");
        printMatrix(b,"B");
        printMatrix(res,"B*A");
        int val = validateProgram2(res);
        if(val==1){
                printf("PROGRAM IS VAlIDATED for B*A\n");
        }else{
                printf("PROGRAM HAS ERRORS!!!\n");
        }
}

if(!strcmp(argv[1],"run")){
createMatrices(a,b,res);
if(PAPI_start_counters(Events,NUM_EVENTS) != PAPI_OK) printf("Nao tem contadores papi\n");
if(PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK) printf("Erro ao ler evento\n");
multMatrices(a,b,res);

if(PAPI_stop_counters(values,NUM_EVENTS) != PAPI_OK) printf("Erro ao parar evento\n");
}

printf("\n---Valores PAPI---\n");
printf("FP_OPS: %lld;\nL3_DCM: %lld;\n",values[0],values[1]);

return 1;
}

