#include <stdio.h>
#include <stdlib.h>

#define SIZE 512


void createMatrices(float **a, float **b){
	int i, j;
	for(i = 0; i < SIZE; i++){
		for(j = 0; j < SIZE; j++){
			a[i][j] = rand() / ((float)RAND_MAX + 1);
			b[i][j] = rand() / ((float)RAND_MAX + 1);
		}
	}
}

void multMatrices(float **a, float **b, float **res){
	int i,j,k;
	for(j = 0; j < SIZE; j++){
		for(i = 0; i < SIZE; i++){
			for(k = 0; k < SIZE; k++){
				res[i][j] += mat1[i][k] * mat2[k][j]; 
			}
		}
	}
}

int main(){
float **a, **b, **res;
int i;

a = malloc(SIZE*sizeof(float*));
b = malloc(SIZE*sizeof(float*));
res = malloc(SIZE*sizeof(float*));

for(i = 0; i < SIZE; i++){
	a[i] = malloc(SIZE*sizeof(float));
	b[i] = malloc(SIZE*sizeof(float));
	res[i] = malloc(SIZE*sizeof(float));
}

createMatrices(a,b);

multMatrices(a,b,res);

return 1;
}
