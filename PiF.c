/*Software Development in the UNIX Environment
Sample C Program

Example C Program to Compute PI Using A Monte Carlo Method

Source code:*/

/* Program to compute Pi using Monte Carlo methods */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#define SEED 35791246


int main(int argc, char**){
   double pi, count_local, x, y, z;
   int rank_W, size_W;

   if(argc != 2){
	printf("Argumentos invalidos");
   }
   
   /* initialize random numbers */
   srand(SEED);
   count = 0;
   
   MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank_W);
	MPI_Comm_size(MPI_COMM_WORLD, &size_W);

	long count, niter = atol(argv[1]) / size_W; 
	
	MPI_Bcast(&niter, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	
	count_local = 0;
	for ( i = 0; i < niter; i++) {
      x = (double)rand()/RAND_MAX;
      y = (double)rand()/RAND_MAX;
      z = x*x+y*y;
      if (z<=1) count_local++;
    }
	
	MPI_Reduce(&count_local, &count, 1, MPI_DOUBLE, MPI_SUM,0, MPI_COMM_WORLD);
	
	if(rank_W == 0){
		pi = 4.0*(double)count/(niter*size_W);
		printf("Valor estimado de pi: %g para %g iteracoes.\n", pi, niter);
	}
	
	MPI_Finalize();
	return 0;
}