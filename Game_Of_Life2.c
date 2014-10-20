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
	int **old, **new, **buf;
	float x;
	
	/* allocate arrays */
	ni = NI + 2; /* add 2 for left and right ghost cells */
	nj = NJ + 2;
	
	printf("START\n");
	MPI_Init(&argc, &argv);
	printf("INIT\n");
	
	int pid, n_proc;
	MPI_Status status;
	// numero de procecos
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);	
	// id do processo
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	printf("\n=== pid: %d\n", pid);

	old = malloc(ni*sizeof(int*));
	new = malloc(ni*sizeof(int*));
	
	for(i=0; i<ni; i++){
	old[i] = malloc(nj*sizeof(int));
	new[i] = malloc(nj*sizeof(int));
	}
	
	//O processo 0 cria a grelha
	if(pid == 0){
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
	}
	
	
	// quantos dados vao para cada processo
	int i_offset = NI / n_proc;
	int offsetL = 0, proc;
	
	// Intervalos
	int i_inicio = pid * i_offset + 1;
	int i_fim = (pid + 1) * i_offset;
	
	// Aloca buffer para enviar novos dados para Master
	buf = malloc(i_offset*sizeof(int*));
	for(i = 0; i < i_offset; i++){
		buf[i] = malloc(NJ*sizeof(int));
	}

	/* time steps */
	for(n=0; n<NSTEPS; n++){
		// Apenas o Master vai fazer isto
		if(pid == 0){
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
		}
		
		// Se for o ultimo tem que se somar o resto
		if (pid+1 == n_proc)
		{i_fim+= NI % n_proc;
		}
		
		// Broadcast da grelha
		MPI_Bcast(&old, ni*nj, MPI_INT, 0, MPI_COMM_WORLD);
		
		//Iterar Ciclos
		for(i=i_inicio; i<=i_fim; i++){
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
				if(n == NSTEPS - 1){
					alive_local++;
				}
			break;
			case 2:
				new[i][j] = old[i][j];
				if(n == NSTEPS - 1){
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
		for(i=0; i<=i_fim; i++){
			for(j=0; j<=NJ; j++){
				old[i][j] = new[i][j];
			}
		}
		if(pid==0){
		// Espero para receber a grelha dos outros
			for (proc = 1; proc < n_proc; proc++){
				if((proc+1) < n_proc){
					MPI_Recv(&buf, i_offset*NJ, MPI_INT, proc, 2, MPI_COMM_WORLD, &status);
					for(i = 0; i <= i_offset; i++){
						for(j = 0; j <= NJ; j++){
							old[i*proc][j] = buf[i][j];
						}
					}
				}
				else{
					offsetL = NI % i_offset;
					MPI_Recv(&buf, offsetL*NJ, MPI_INT, proc, 2, MPI_COMM_WORLD, &status);
					for(i = 0; i <= offsetL; i++){
						for(j = 0; j <= NJ; j++){
							old[i*proc][j] = buf[i][j];
						}
					}
				}
			}
		}
		else // Sou Slave e envio a minha nova grelha
		{
			for(i = 0; i <= i_offset; i++){
				for(j = 0; j <= NJ; j++){
					buf[i][j] = new[i+i_inicio][j];
				}
			}
			MPI_Send(&buf, i_offset*NJ, MPI_INT, 0, 2, MPI_COMM_WORLD);
		}
		
		if(n == NSTEPS - 1){
			MPI_Reduce(&alive_local, &alive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
		}
	}
	
	if(pid == 0){
		printf("\nNumber of live cells = %d\n", alive);
	}
	
	free(old);
	free(new);
	free(buf);
	
	MPI_Finalize();
	
	return 0;
	}