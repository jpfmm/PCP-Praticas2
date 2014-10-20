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
	
	printf("START\n");
	MPI_Init(&argc, &argv);
	printf("INIT\n");
	
	int i, j, n, im, ip, jm, jp, ni, nj, nsum, isum;
	int **old, **new,buf;
	float x;
	
	/* allocate arrays */
	ni = NI + 2; /* add 2 for left and right ghost cells */
	nj = NJ + 2;
	
	int pid, n_proc;
	MPI_Status status;
	// numero de procecos
	MPI_Comm_size(MPI_COMM_WORLD, &n_proc);	
	// id do processo
	MPI_Comm_rank(MPI_COMM_WORLD, &pid);
	printf("\n=== pid: %d\n", pid);
	
	// Gera no Master
	if(pid==0){
	old = malloc(ni*sizeof(int*));
	new = malloc(ni*sizeof(int*));
	//temp = malloc(ni*sizeof(int*));
	
	for(i=0; i<ni; i++){
	old[i] = malloc(nj*sizeof(int));
	new[i] = malloc(nj*sizeof(int));
	//temp[i] = malloc(nj*sizeof(int));
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
	printf("FIM ALOC\n");
	}
	
	// quantos dados vao para cada processo
	int i_offset = NI / n_proc;
	int j_offset = NJ / n_proc;
	int off = 0, proc;
	
	// Intervalos
	int i_inicio = pid * i_offset + 1;
	int i_fim = (pid + 1) * i_offset;
	
	// Enviar o endereco do old e new
	if(pid == 0){
	for (proc = 1; proc < n_proc; proc++){
		MPI_Send(&old, 1, MPI_INT, proc, 10, MPI_COMM_WORLD);
		MPI_Send(&new, 1, MPI_INT, proc, 11, MPI_COMM_WORLD);}
	}else{
		MPI_Recv(&old, 1, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
		MPI_Recv(&new, 1, MPI_INT, 0, 11, MPI_COMM_WORLD, &status);
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
		// Enviar a todos que podem prosseguir
		for (proc = 1; proc < n_proc; proc++)
		MPI_Send(&buf, 1, MPI_INT, proc, 1, MPI_COMM_WORLD);
		}
	else // Sou Slave e espero pelo Master para poder prosseguir
	{
	MPI_Recv(&buf, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
	}
	// Se for o ultimo tem que se somar o resto
	if (pid+1 == n_proc)
	{i_fim+= NI % n_proc;
	}
	
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
		break;
		case 2:
			new[i][j] = old[i][j];
		break;
		default:
			new[i][j] = 0;
		}
	}
	}
	// Feito o calculo vou gravar os dados
	// Copio os meus valores
	for(i=0; i<=i_fim; i++)
	for(j=0; j<=NJ; j++){
		old[i][j] = new[i][j];
	}
	if(pid==0){
	// Espero pelos outros
	for (proc = 1; proc < n_proc; proc++)
	MPI_Recv(&buf, 1, MPI_INT, proc, 2, MPI_COMM_WORLD, &status);
	}
	else // Sou Slave e espero pelo Master para poder prosseguir
	{
	MPI_Send(&buf, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	}
	}
	/*
	// Receber as coisas
	MPI_Send(&bin_count, 1, MPI_INT, proc, 1, MPI_COMM_WORLD);
	*/
	MPI_Finalize();
	/* Iterations are done; sum the number of live cells */
	isum = 0;
	for(i=1; i<=NI; i++){
	for(j=1; j<=NJ; j++){
	isum = isum + new[i][j];
	}
	}
	printf("\nNumber of live cells = %d\n", isum);
	return 0;
	}