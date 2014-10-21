/***********************
Conway Game of Life
serial version
************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> //MPI
#define NI 200 /* array sizes */
#define NJ 200
#define NSTEPS 500 /* number of time steps */

int main(int argc, char *argv[]) {
	
	int i, j, n, im, ip, jm, jp, ni, nj, nsum, alive = 0, alive_local = 0;
	int **old, **new;
	float x;
	
	MPI_Init(&argc, &argv);
	
	/* allocate arrays */
	ni = NI + 2; /* add 2 for left and right ghost cells */
	nj = NJ + 2;
	
	int pid, n_proc;
	MPI_Status status;
	// numero de procecos
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);	
	// id do processo
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);

	old = malloc(ni*sizeof(int*));
	new = malloc(ni*sizeof(int*));
	
	for(i=0; i<ni; i++){
	old[i] = malloc(nj*sizeof(int));
	new[i] = malloc(nj*sizeof(int));
	}
	
	
	/* seed */
	srand(time(0));

	/* initialize elements of old to 0 or 1 */
	for(i=1; i<=NI; i++){
	for(j=1; j<=NJ; j++){
		x = rand()/((float)RAND_MAX + 1);
		if(x<0.5){
			old[i][j] = 0;
		} else {
			old[i][j] = 1;
		}
	}
	}
	
	//Quantas simulacoes cada processo faz
	int carga = NSTEPS / n_proc;
	
	// Se for o ultimo faz o resto das iteracoes
		if (pid+1 == n_proc)
		{carga+= NSTEPS % n_proc;
		}
	
	// Efectua os calculos respetivos
	for(n=0; n<carga; n++){
		
		/* corner boundary conditions */
		old[0][0] = old[NI][NJ];
		old[0][NJ+1] = old[NI][1];
		old[NI+1][NJ+1] = old[1][1];
		old[NI+1][0] = old[1][NJ];
		/* left-right boundary conditions */
		for(i=1; i<=NI; i++){
		old[i][0] = old[i][NJ];
		old[i][NJ+1] = old[i][1];
		}
		/* top-bottom boundary conditions */
		for(j=1; j<=NJ; j++){
		old[0][j] = old[NI][j];
		old[NI+1][j] = old[1][j];
		}
		
		//Iterar Ciclos
		for(i=1; i<=NI; i++){
		for(j=0; j<=NJ; j++){
			im = i-1;
			ip = i+1;
			jm = j-1;
			jp = j+1;
			nsum = old[im][jp] + old[i][jp] + old[ip][jp]
			+ old[im][j ] + old[ip][j ]
			+ old[im][jm] + old[i][jm] + old[ip][jm];
			switch(nsum){
			case 3:
				new[i][j] = 1;
				if(n == (carga - 1)){
					alive_local++;
				}
			break;
			case 2:
				new[i][j] = old[i][j];
				if(n == (carga - 1)){
					alive_local += old[i][j];
				}
			break;
			default:
				new[i][j] = 0;
			}
		}
		}
		// Feito o calculo vou gravar os dados
		// Copio os meus valores
		for(i=0; i<=NI; i++){
			for(j=0; j<=NJ; j++){
				old[i][j] = new[i][j];
			}
		}
	}
	
	MPI_Reduce(&alive_local, &alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	printf("Hi, PID: %d here, so i got %d guys alive\n", alive_local);
	
	if(pid == 0){
		printf("\nNumber of live cells = %d, for a %dx%d grid with %d iterations\n", alive, NI, NJ, NSTEPS);
	}
	
	free(old);
	free(new);
	
	MPI_Finalize();
	
	return 0;
	}