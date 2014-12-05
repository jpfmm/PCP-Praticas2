#include <stdio.h>
#include <stdlib.h>

#define SIZE 512


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

int main(){
float a[SIZE][SIZE], b[SIZE][SIZE], res[SIZE][SIZE];
int i;

createMatrices(a,b);

multMatrices(a,b,res);

return 1;
}
