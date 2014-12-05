#include <stdio.h>
#include <stdlib.h>

#define SIZE 512


void createValidMatrices(float a[SIZE][SIZE], float b[SIZE][SIZE]){
	int i,j;
	for(i = 0; i < SIZE; i++){
		for(j = 0; j < SIZE; j++){
			a[i][j] = rand() / ((float)RAND_MAX + 1);
			b[i][j] = 1;
		}
	}
}

void createMatrices(float a[SIZE][SIZE], float b[SIZE][SIZE]){
	int i, j;
	for(i = 0; i < SIZE; i++){
		for(j = 0; j < SIZE; j++){
			a[i][j] = rand() / ((float)RAND_MAX + 1);
			b[i][j] = rand() / ((float)RAND_MAX + 1);
		}
	}
}

void multMatrices(float a[SIZE][SIZE], float b[SIZE][SIZE], float res[SIZE][SIZE]){
	int i,j,k;
	for(j = 0; j < SIZE; j++){
		for(i = 0; i < SIZE; i++){
			for(k = 0; k < SIZE; k++){
				res[i][j] += a[i][k] * b[k][j]; 
			}
		}
	}
}

int validadeProgram1(float res[SIZE][SIZE]){
	int i,j, result=1;
	float ref=0, aux; 
	for(i = 0; i < SIZE; i++){
		aux = 0;
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

int validadeProgram2(float res[SIZE][SIZE]){
	int i,j, result=1;
	float ref=0, aux; 
	for(j = 0; j < SIZE; j++){
		aux = 0;
		for(i = 0; i < SIZE; i++){
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

int main(int argc, char *argv[]){
float a[SIZE][SIZE], b[SIZE][SIZE], res[SIZE][SIZE];
int i;

if(argc!=2){
	printf("Argumentos inválidos\n USAGE: prog_name FLAG\n");
	return 1;
}

if(argv[1]=="VALIDATE1"){
	createValidMatrices(a,b);
	multMatrices(a,b,res);
	int val = validateProgram1(a,res);
	if(val==1){
		printf("PROGRAM IS VAlIDATED for A*B\n");
	}else{
		preintf("PROGRAM HAS ERRORS!!!");
	}
}

if(argv[1]=="VALIDATE2"){
	createValidMatrices(a,b);
	multMatrices(b,a,res);
	int val = validateProgram2(a,res);
	if(val==1){
		printf("PROGRAM IS VAlIDATED for B*A\n");
	}else{
		preintf("PROGRAM HAS ERRORS!!!");
	}
}

if(argv[1]=="RUN"){
createMatrices(a,b);

multMatrices(a,b,res);
}
return 1;
}
