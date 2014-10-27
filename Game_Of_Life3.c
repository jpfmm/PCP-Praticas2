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
	int **new, **old, *buf_send, *buf_recv;
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
	
	MPI_Request send1, send2;
	MPI_Request recieve1, recieve2;
	MPI_Status status;
	
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
	
	buf_send = malloc(nj*sizeof(int));
 	buf_recv = malloc(nj*sizeof(int));
	
	/* seed */
	// srand(time(0));

	/* initialize elements of old to 0 or 1 */
	for(i=1; i<(offset-1); i++){
	for(j=1; j<=NJ; j++){
		x = rand()/((float)RAND_MAX + 1);
		if(x<0.5){
			old[i][j] = 0;
		} else {
			old[i][j] = 1;
		}
	}
	}

	
	//Inicializa as comunicações
	if(pid == 0){
		MPI_Send_init(&buf_send, nj, MPI_INT, n_proc-1, 2, MPI_COMM_WORLD, &send1);
	}else{
		MPI_Send_init(&buf_send, nj, MPI_INT, pid-1, 2, MPI_COMM_WORLD, &send1);
	}
	if(pid == (n_proc-1)){
		MPI_Send_init(&buf_send, nj, MPI_INT, 0, 3, MPI_COMM_WORLD, &send2);
	}else{
		MPI_Send_init(&buf_send, nj, MPI_INT, pid+1, 3, MPI_COMM_WORLD, &send2);
	}
	
	if(pid == (n_proc-1)){
		MPI_Recv_init(&buf_recv, nj, MPI_INT, 0, 2, MPI_COMM_WORLD, &recieve1);
	}else{
		MPI_Recv_init(&buf_recv, nj, MPI_INT, pid+1, 2, MPI_COMM_WORLD, &recieve1);
	}
	if(pid == 0){
		MPI_Recv_init(&buf_recv, nj, MPI_INT, n_proc-1, 3, MPI_COMM_WORLD, &recieve2);
	}else{
		MPI_Recv_init(&buf_recv, nj, MPI_INT, pid-1, 3, MPI_COMM_WORLD, &recieve2);
	}
	
	
	// Efectua os calculos respetivos
	for(n = 0; n < NSTEPS; n++){
		
		//Actualizo os ghosts da direita e da esquerda
		for(i=1; i<offset-1; i++){
			old[i][0] = old[i][NJ];
			old[i][NJ+1] = old[i][1];
		}
		
		//Agora so tenho que actualizar a minha grelha e enviar a primeira e a ultima linha
		//para o processo anterior e para o posterior, respetivamente
		for(i=0; i<nj; i++){
			buf_send[i] = old[0][i];
		}
		MPI_Start(&send1);
		for(i=0; i<nj; i++){
			buf_send[i] = old[ni-1][i];
		}
		MPI_Start(&send2);
		
		//Tenho que receber e gravar
		MPI_Start(&recieve1);
		for(i = 0; i < nj; i++){
			old[offset-1][i] = buf_recv[i];
		}
		
		MPI_Start(&recieve2);
		for(i = 0; i < nj; i++){
			old[0][i] = buf[i];
		}
		//Garanto que os sends e recieves estao feitos
		MPI_Wait (&send1, status);
		MPI_Wait (&send2, status);
		MPI_Wait (&recieve1, status);
		MPI_Wait (&recieve2, status);
		
		//Iterar Ciclos
		for(i=1; i<offset-1; i++){
		for(j=1; j<=NJ; j++){
			im = i-1;
			ip = i+1;
			jm = j-1;
			jp = j+1;
			nsum = old[im][jp] + old[i][jp] + old[ip][jp]
			+ old[im][j ] 				+ old[ip][j ]
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
	}
	
	MPI_Reduce(&alive_local, &alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	
	if(pid == 0){
		printf("\nNumber of live cells = %d, for a %dx%d grid with %d iterations\n", alive, NI, NJ, NSTEPS);
	}
	
	//Eliminar as conexoes
	MPI_Request_free (&send1);
	MPI_Request_free (&send2);
	MPI_Request_free (&recieve1);
	MPI_Request_free (&recieve2);
	
	//Eliminar alocacoes
	free(old);
	free(new);
	free(buf_send);
	free(buf_recv);
	
	MPI_Finalize();
	
	return 0;
	}