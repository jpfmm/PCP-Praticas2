/***********************
Conway Game of Life
serial version

Autor: João Magalhães pg27762

Compile:   mpicc -o GameOfLife Game_Of_Life.c
Run:       ./GameOfLife 
************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h> //MPI
#define NI 200 /* array sizes */
#define NJ 200
#define NSTEPS 500 /* number of time steps */

int main(int argc, char *argv[]) {
	
	int i, j, n, im, ip, jm, jp, ni, nj, nsum, alive = 0, alive_local = 0;
	int **ori, **new, **old, *buf;
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
	
	int offset = NI / n_proc + 2;

	if (pid+1 == n_proc){
				offset+= NI % n_proc;
			}
			
	old = malloc(offset*sizeof(int*));
	new = malloc(offset*sizeof(int*));
	
	for(i=0; i<offset; i++){
	old[i] = malloc(nj*sizeof(int));
	new[i] = malloc(nj*sizeof(int));
	}
	
	buf = malloc(nj*sizeof(int));
 	
	if(pid == 0){
		ori = malloc(ni*sizeof(int*));
		
		for(i=0; i<ni; i++){
		ori[i] = malloc(nj*sizeof(int));
		}
	
		/* seed */
		srand(time(0));

		/* initialize elements of old to 0 or 1 */
		for(i=1; i<=NI; i++){
		for(j=1; j<=NJ; j++){
			x = rand()/((float)RAND_MAX + 1);
			if(x<0.5){
				ori[i][j] = 0;
			} else {
				ori[i][j] = 1;
			}
		}
		}
		
		//Calcula os cantos
		/* corner boundary conditions */
		ori[0][0] = ori[NI][NJ];
		ori[0][NJ+1] = ori[NI][1];
		ori[NI+1][NJ+1] = ori[1][1];
		ori[NI+1][0] = ori[1][NJ];
		/* left-right boundary conditions */
		for(i=1; i<=NI; i++){
		ori[i][0] = ori[i][NJ];
		ori[i][NJ+1] = ori[i][1];
		}
		/* top-bottom boundary conditions */
		for(j=1; j<=NJ; j++){
		ori[0][j] = ori[NI][j];
		ori[NI+1][j] = ori[1][j];
		}
		
		//passa para a sua grelha de calculo
		for(i = 0; i < offset; i++){
			for(j = 0; j < nj; j++){
				old[i][j] = ori[i][j];
			}
		}
		
		for(i = 1; i < n_proc - 1; i++){
			MPI_Send(&ori + (i*offset)*nj, offset, MPI_INT, i, 1, MPI_COMM_WORLD);
		}
		
		MPI_Send(&ori + (i*offset)*nj, (offset + (NI % n_proc))*nj, MPI_INT, i, 1, MPI_COMM_WORLD);
		
		free(ori);
	}
	else{
		MPI_Recv(&old, offset*nj, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	
	
	// Efectua os calculos respetivos
	for(n = 0; n < NSTEPS; n++){
		
		//Iterar Ciclos
		for(i=1; i<offset-1; i++){
		for(j=1; j<=NJ; j++){
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
				if(n == (NSTEPS - 1)){
					alive_local++;
				}
			break;
			case 2:
				new[i][j] = old[i][j];
				if(n == (NSTEPS - 1)){
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
		for(i=1; i<offset-1; i++){
			for(j=0; j<=NJ; j++){
				old[i][j] = new[i][j];
			}
		}
		
		//Actualizo os ghosts da direita e da esquerda
		for(i=1; i<offset-1; i++){
			old[i][0] = old[i][NJ];
			old[i][NJ+1] = old[i][1];
		}
		
		//Agora so tenho que actualizar a minha grelha e enviar a primeira e a ultima linha
		//para o processo anterior e para o posterior, respetivamente
		
		MPI_Send(&old + offset, nj, MPI_INT, pid-1, 2, MPI_COMM_WORLD);
		MPI_Send(&old + offset*nj, nj, MPI_INT, pid+1, 3, MPI_COMM_WORLD);
		
		//Tenho que receber e gravar
		
		MPI_Recv(&buf, nj, MPI_INT, pid+1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(i = 0; i < nj; i++){
			old[offset-1][i] = buf[i];
		}
		MPI_Recv(&buf, nj, MPI_INT, pid-1, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for(i = 0; i < nj; i++){
			old[0][i] = buf[i];
		}
		
		
	}
	
	MPI_Reduce(&alive_local, &alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	
	if(pid == 0){
		printf("\nNumber of live cells = %d, for a %dx%d grid with %d iterations\n", alive, NI, NJ, NSTEPS);
	}
	
	free(old);
	free(new);
	free(buf);
	
	MPI_Finalize();
	
	return 0;
	}